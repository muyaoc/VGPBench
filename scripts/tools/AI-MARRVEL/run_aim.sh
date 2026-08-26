#!/usr/bin/env bash
set -euo pipefail

# Record script start time
SCRIPT_START_TIME=$(date +%s)

############################################
# Config
############################################
# Script directory (this file lives in scripts/tools/AI-MARRVEL/),
# so main.nf and nextflow.config sit right next to it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Base directory of the benchmark + tool data.
# Set GENETOOLS_BASE before running, e.g.:
#   export GENETOOLS_BASE=/path/to/genetools
if [ -z "${GENETOOLS_BASE:-}" ]; then
    echo "ERROR: set GENETOOLS_BASE to your project root (e.g. export GENETOOLS_BASE=/path/to/genetools)" >&2
    exit 1
fi

AIM_REF_DIR="$GENETOOLS_BASE/aim/aim_data"
VCF_ROOT="$GENETOOLS_BASE/genotype"
HPO_ROOT="$GENETOOLS_BASE/phenotype"

OUT_ROOT="$GENETOOLS_BASE/aim/results"
STORE_ROOT="$GENETOOLS_BASE/aim/store"

PIPELINE_MAIN="${SCRIPT_DIR}/main.nf"
PROFILE="docker"

JOBS=${JOBS:-4}
MAX_CPUS=${MAX_CPUS:-2}
MAX_MEM=${MAX_MEM:-"20.GB"}
############################################

export NXF_DISABLE_CHECK_UPDATES=true
export NXF_OFFLINE=true
# 限制 Nextflow 本身的资源使用
export NXF_JVM_OPTS="-Xms1g -Xmx4g"

log(){ echo "[$(date '+%F %T')] $*"; }

# Log initial configuration
log "=============================================="
log "CONFIGURATION"
log "=============================================="
log "Parallel jobs: $JOBS"
log "CPUs per job: $MAX_CPUS"
log "Memory per job: $MAX_MEM"
log "Total CPUs: $((JOBS * MAX_CPUS))"
log "Host: $(hostname)"
log "Pipeline: $PIPELINE_MAIN"
log "Profile: $PROFILE"
log "=============================================="

# run_id from vcf filename: patient_<run_id>.vcf(.gz)
runid_from_vcf() {
  local bn
  bn="$(basename "$1")"
  bn="${bn%.vcf.gz}"
  bn="${bn%.vcf}"
  bn="${bn#patient_}"
  echo "$bn"
}

# outdir name: take "caseNNN" from run_id like "case11_HG00525" -> "case11"
outname_from_runid() {
  local rid="$1"
  if [[ "$rid" =~ ^(case[0-9]+) ]]; then
    echo "${BASH_REMATCH[1]}"
  else
    # fallback: use full run_id
    echo "$rid"
  fi
}

# hpo path from run_id: <HPO_ROOT>/<GRChXX>/<run_id>_hpos.txt
hpo_path() {
  local refdir="$1" rid="$2"
  echo "$HPO_ROOT/$refdir/${rid}_hpos.txt"
}

# Generate possible final expanded paths (to cover all current and previous nested structures)
get_candidate_paths() {
  local rid="$1"
  local shortname
  shortname="$(outname_from_runid "$rid")"

  # Current structure (1 layer using run_id) e.g. case5_HG00448
  echo "$OUT_ROOT/$rid/prediction/final_matrix_expanded/${rid}.expanded.csv.gz"
  # Previous structure (2 layers using shortname/run_id) e.g. case1/case1_HG00407
  echo "$OUT_ROOT/$shortname/$rid/prediction/final_matrix_expanded/${rid}.expanded.csv.gz"
  # Previous structure (2 layers using run_id/run_id) e.g. case1006_HG03653/case1006_HG03653
  echo "$OUT_ROOT/$rid/$rid/prediction/final_matrix_expanded/${rid}.expanded.csv.gz"
  # Extra fallback (1 layer using shortname)
  echo "$OUT_ROOT/$shortname/prediction/final_matrix_expanded/${rid}.expanded.csv.gz"
}

# Output path for printing missing/bad warnings matching the current logic
expected_output_path() {
  local rid="$1"
  echo "$OUT_ROOT/$rid/prediction/final_matrix_expanded/${rid}.expanded.csv.gz"
}

# Check if ANY of the candidate paths are completed correctly
is_done() {
  local rid="$1"
  local paths
  mapfile -t paths < <(get_candidate_paths "$rid")
  
  for f in "${paths[@]}"; do
    if [[ -f "$f" && -s "$f" ]]; then
      if gzip -t "$f" >/dev/null 2>&1; then
        return 0
      fi
    fi
  done
  return 1
}

############################################
# sanity
############################################
command -v nextflow >/dev/null 2>&1 || { log "ERROR: nextflow not found"; exit 1; }
command -v gzip >/dev/null 2>&1 || { log "ERROR: gzip not found"; exit 1; }
[[ -f "$PIPELINE_MAIN" ]] || { log "ERROR: main.nf not found: $PIPELINE_MAIN"; exit 1; }

mkdir -p "$OUT_ROOT/_work" "$OUT_ROOT/_logs" "$STORE_ROOT"

############################################
# collect vcf list
############################################
mapfile -t VCF_LIST < <(
  find "$VCF_ROOT/GRCh37" -maxdepth 1 -type f -name "patient_*.vcf.gz" -o -name "patient_*.vcf" 2>/dev/null
  find "$VCF_ROOT/GRCh38" -maxdepth 1 -type f -name "patient_*.vcf.gz" -o -name "patient_*.vcf" 2>/dev/null
)

# prefer .vcf.gz when both exist
declare -A VCF_BY_RUNID REF_BY_RUNID
for v in "${VCF_LIST[@]}"; do
  [[ -f "$v" ]] || continue
  rid="$(runid_from_vcf "$v")"
  refdir="GRCh37"
  refver="hg19"
  [[ "$v" == *"/GRCh38/"* ]] && refdir="GRCh38" && refver="hg38"

  prev="${VCF_BY_RUNID[$rid]:-}"
  if [[ -z "$prev" ]]; then
    VCF_BY_RUNID["$rid"]="$v"
    REF_BY_RUNID["$rid"]="$refver|$refdir"
  else
    # prefer gz
    if [[ "$prev" =~ \.vcf$ ]] && [[ "$v" =~ \.vcf\.gz$ ]]; then
      VCF_BY_RUNID["$rid"]="$v"
      REF_BY_RUNID["$rid"]="$refver|$refdir"
    fi
  fi
done

mapfile -t RUNIDS < <(printf '%s\n' "${!VCF_BY_RUNID[@]}" | sort)

############################################
# run with simple concurrency
############################################
running=0
skipped=0
launched=0

run_one() {
  local rid="$1"
  local vcf="$2"
  local refver="$3"
  local refdir="$4"

  local outname outdir storedir workdir logfile hpo
  outname="$rid"
  outdir="$OUT_ROOT/$outname"
  storedir="$STORE_ROOT/$outname"
  workdir="$OUT_ROOT/_work/$rid"
  logfile="$OUT_ROOT/_logs/${rid}.log"
  hpo="$(hpo_path "$refdir" "$rid")"

  mkdir -p "$outdir" "$storedir" "$workdir"

  # skip if done
  if is_done "$rid"; then
    log "SKIP $rid (done)"
    return 0
  fi

  # require hpo
  if [[ ! -f "$hpo" ]]; then
    log "SKIP $rid (missing HPO: $hpo)"
    return 0
  fi

  # Log start time and run
  local start_time end_time duration
  start_time=$(date +%s)
  log "START $rid  ref=$refver  cpus=$MAX_CPUS  mem=$MAX_MEM"

  nextflow run "$PIPELINE_MAIN" -profile "$PROFILE" -work-dir "$workdir" \
    -process.maxCpus="$MAX_CPUS" -process.maxMemory="$MAX_MEM" \
    --ref_dir "$AIM_REF_DIR" \
    --input_vcf "$vcf" \
    --input_hpo "$hpo" \
    --ref_ver "$refver" \
    --run_id "$rid" \
    --outdir "$outdir" \
    --storedir "$storedir" \
    >"$logfile" 2>&1
  local rc=$?

  # Calculate and log duration
  end_time=$(date +%s)
  duration=$((end_time - start_time))
  local duration_min=$((duration / 60))
  local duration_sec=$((duration % 60))

  # if rc=0 but final missing/bad, mark in log
  if ! is_done "$rid"; then
    log "FAIL  $rid  duration=${duration_min}m${duration_sec}s  rc=$rc  (output missing)"
    echo "[WARN] finished but final expanded.gz missing/bad: $(expected_output_path "$rid")" >> "$logfile"
    return 2
  fi

  log "DONE  $rid  duration=${duration_min}m${duration_sec}s  rc=$rc"
}

for rid in "${RUNIDS[@]}"; do
  vcf="${VCF_BY_RUNID[$rid]}"
  IFS='|' read -r refver refdir <<< "${REF_BY_RUNID[$rid]}"

  # 先判断是否 skip（避免占用并发槽）
  if is_done "$rid"; then
    skipped=$((skipped+1))
    log "SKIP $rid (done)"
    continue
  fi

  # 后台跑
  (
    set +e
    run_one "$rid" "$vcf" "$refver" "$refdir"
    exit $?
  ) &

  launched=$((launched+1))
  running=$((running+1))

  if (( running >= JOBS )); then
    set +e
    wait -n
    set -e
    running=$((running-1))
  fi
done

while (( running > 0 )); do
  set +e
  wait -n
  set -e
  running=$((running-1))
done

# Final summary
SCRIPT_END_TIME=$(date +%s)
SCRIPT_DURATION=$((SCRIPT_END_TIME - SCRIPT_START_TIME))
SCRIPT_DURATION_MIN=$((SCRIPT_DURATION / 60))
SCRIPT_DURATION_SEC=$((SCRIPT_DURATION % 60))

log "=============================================="
log "SUMMARY"
log "=============================================="
log "Total launched: $launched"
log "Total skipped (done): $skipped"
log "Concurrency: $JOBS parallel jobs"
log "Per-job config: $MAX_CPUS CPUs, $MAX_MEM memory"
log "Host: $(hostname)"
log "Total runtime: ${SCRIPT_DURATION_MIN}m${SCRIPT_DURATION_SEC}s"
log "Logs: $OUT_ROOT/_logs"
log "=============================================="