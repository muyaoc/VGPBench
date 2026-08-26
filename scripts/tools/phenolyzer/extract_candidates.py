#!/usr/bin/env python3
"""
extract_candidates.py

Extract candidate gene lists from raw joint-called VCFs using ANNOVAR.

Strategy (per VCF):
  1. Decompress .vcf.gz with Python gzip  (no gunzip dependency)
  2. Parse INFO field → population AF; find rare loci (AF < threshold)
     AND where the sample has non-ref genotype (GT != 0/0)
  3. Write a pre-filtered VCF containing ONLY those rare loci
  4. convert2annovar.pl  → avinput  (on filtered VCF — fast)
  5. table_annovar.pl    → refGene  (on filtered avinput — fast)
  6. Parse multianno.txt → extract gene symbols
  7. Write <caseID>_<SampleID>_genes.txt

Output: <OUT_DIR>/<caseID>_<SampleID>_genes.txt  (one gene symbol per line)

Requirements:
  - ANNOVAR (on PATH; e.g. convert2annovar.pl, table_annovar.pl)
      with humandb/hg19_refGene.txt
  - Perl (with BioPerl, for ANNOVAR)
  - Python 3.8+ with built-in gzip module

Usage:
    python extract_candidates.py                          # all cases
    python extract_candidates.py --workers 4            # parallel
    python extract_candidates.py --dry-run               # preview
    python extract_candidates.py --cases case1 case2   # specific cases
"""

import os
import re
import sys
import time
import logging
import argparse
import subprocess
import multiprocessing
from pathlib import Path
from typing import Optional, List, Tuple, Dict, Set

# ── Logging ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)-7s %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("extract_candidates")

# ── Paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent.resolve()
# Perl binary: set PERL_BIN if perl is not on PATH (e.g. on Windows:
#   set PERL_BIN=D:\Strawberry\perl\bin\perl.exe)
STRAWBERRY_PERL = Path(os.environ.get("PERL_BIN", "perl"))
ANNOVAR_DIR = SCRIPT_DIR / "annovar"
ANNOVAR_DB = ANNOVAR_DIR / "humandb"
VCF_DIR = SCRIPT_DIR / "genotype" / "GRCh37"
OUT_DIR = SCRIPT_DIR / "phenolyzer_candidates"

# Defaults (overridden by CLI args)
DEFAULT_AF_THRESHOLD = 0.01
DEFAULT_MIN_GQ = 10

# ── gzip lazy import ──────────────────────────────────────────────────────────
_gzip_mod: Optional[object] = None


def _ensure_gzip():
    global _gzip_mod
    if _gzip_mod is None:
        import gzip
        _gzip_mod = gzip


# ── Core functions ────────────────────────────────────────────────────────────

def _parse_vcf_find_rare_loci(
    vcf_path: Path,
    rare_af: float,
    min_gq: int,
) -> Tuple[Dict[str, float], Set[str]]:
    """
    Read a decompressed VCF and build:
      af_lookup : {locus_key → float(AF)}  from INFO field
      rare_loci : set of locus keys where AF < threshold AND GT != 0/0

    Locus key format: "chr1" (ANNOVAR style)
    """
    af_lookup: Dict[str, float] = {}
    rare_loci: Set[str] = set()

    opener = _gzip_mod.open if vcf_path.suffix == ".gz" else open
    fmt_idx = -1
    sample_idx = -1

    with opener(vcf_path, "rt", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip()
            if line.startswith("#CHROM"):
                cols = line.split("\t")
                for i, c in enumerate(cols):
                    if c == "FORMAT":
                        fmt_idx = i
                    elif i > fmt_idx and fmt_idx >= 0:
                        sample_idx = i
                        break
                continue
            if line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 8:
                continue

            chrom = fields[0].strip()
            pos = fields[1].strip()
            info_str = fields[7]

            # Parse AF from INFO field
            af = None
            for tok in info_str.split(";"):
                if tok.startswith("AF="):
                    try:
                        af = float(tok[3:])
                        break
                    except ValueError:
                        pass

            # Normalise chromosome key (ANNOVAR uses chr1, not 1)
            key_chr = f"chr{chrom}" if not chrom.startswith("chr") else chrom
            key_num = chrom
            key_full_chr = f"{key_chr}:{pos}"
            key_full_num = f"{key_num}:{pos}"

            # Skip common variants (AF >= threshold)
            # NOTE: variants without AF field (e.g. HGVS-only patient-specific variants)
            # are kept and treated as rare (AF=0) — they will be filtered by ANNOVAR later
            if af is not None and af >= rare_af:
                continue

            if af is not None:
                af_lookup[key_full_chr] = af
                af_lookup[key_full_num] = af

            # Parse per-sample genotype
            if sample_idx >= 0 and sample_idx < len(fields):
                fmt = fields[fmt_idx].split(":") if fmt_idx >= 0 else []
                sample_data = fields[sample_idx].split(":")
                if not fmt or not sample_data:
                    continue

                gt_idx = next((j for j, f in enumerate(fmt) if f == "GT"), None)
                gq_idx = next((j for j, f in enumerate(fmt) if f == "GQ"), None)

                has_alt = False
                if gt_idx is not None and gt_idx < len(sample_data):
                    gt = sample_data[gt_idx]
                    gt_clean = gt.replace("|", "/")
                    # 0/0, ./., .  → reference
                    # 0/1, 1/1, 1|0, etc. → alternate
                    if gt_clean not in ("0/0", "./.", ".") and not gt_clean.startswith("0/0"):
                        has_alt = True

                if has_alt and gq_idx is not None and gq_idx < len(sample_data):
                    try:
                        gq = int(sample_data[gq_idx])
                        if gq < min_gq:
                            has_alt = False
                    except ValueError:
                        pass

                if has_alt:
                    rare_loci.add(f"{key_chr}:{pos}")
                    rare_loci.add(f"{key_num}:{pos}")

    return af_lookup, rare_loci


def _write_filtered_vcf(
    vcf_path: Path,
    rare_loci: Set[str],
    out_path: Path,
) -> int:
    """
    Write a VCF containing only the rare loci.
    Returns number of variant lines written.
    """
    opener = _gzip_mod.open if vcf_path.suffix == ".gz" else open
    n_written = 0

    with opener(vcf_path, "rt", encoding="utf-8", errors="replace") as fin:
        with open(out_path, "w", encoding="utf-8") as fout:
            for raw in fin:
                if raw.startswith("#"):
                    fout.write(raw)
                    continue
                fields = raw.rstrip().split("\t")
                if len(fields) < 2:
                    continue
                chrom = fields[0].strip()
                pos = fields[1].strip()
                key_chr = f"chr{chrom}" if not chrom.startswith("chr") else chrom
                key_num = chrom
                if f"{key_chr}:{pos}" in rare_loci or f"{key_num}:{pos}" in rare_loci:
                    fout.write(raw)
                    n_written += 1

    return n_written


def parse_case_id_from_vcf(vcf_path: Path) -> Tuple[str, str]:
    """patient_case1000_HG03643.vcf.gz → ('case1000', 'HG03643')"""
    # Strip .gz and .vcf suffixes first
    name = vcf_path.name
    name = re.sub(r'\.vcf\.gz$', '', name, flags=re.IGNORECASE)
    name = re.sub(r"^patient_", "", name)
    parts = name.rsplit("_", 1)
    return (parts[0], parts[1]) if len(parts) == 2 else (name, "")


def extract_genes_from_multianno(
    multianno_path: Path,
    rare_loci: Set[str],
    af_lookup: Dict[str, float],
    rare_af: float,
) -> List[str]:
    """
    Parse ANNOVAR .hg19_multianno.txt and extract gene symbols for:
      (a) rare loci (AF < threshold) AND non-ref genotype, OR
      (b) coding/splicing impact regardless of AF
    """
    genes: Dict[str, bool] = {}
    gene_col = -1
    func_col = -1
    chr_col = -1
    pos_col = -1

    try:
        with open(multianno_path, encoding="utf-8", errors="replace") as fh:
            raw = fh.readline().rstrip()
            raw = re.sub(r"^(NULL|\xef\xbf\xbd)\t?", "", raw)
            cols = [c.strip().strip('"') for c in raw.split("\t")]
            for i, c in enumerate(cols):
                if c == "Gene.refGene":
                    gene_col = i
                elif c == "Func.refGene":
                    func_col = i
                elif c.lower() == "chr":
                    chr_col = i
                elif c.lower() == "start":
                    pos_col = i

            if gene_col < 0:
                log.warning("Gene.refGene not found in %s", multianno_path)
                return []

        with open(multianno_path, encoding="utf-8", errors="replace") as fh:
            fh.readline()  # skip header
            for line in fh:
                if not line.strip():
                    continue
                line = re.sub(r"^NULL\t?", "", line)
                cols = [c.strip().strip('"') for c in line.rstrip().split("\t")]

                if gene_col < 0 or gene_col >= len(cols):
                    continue

                gene_str = cols[gene_col]
                if not gene_str or gene_str in (".", "-"):
                    continue
                func_str = cols[func_col] if func_col >= 0 and func_col < len(cols) else ""

                # Build locus key from Chr + Start
                if chr_col >= 0 and pos_col >= 0:
                    chr_v = cols[chr_col].strip()
                    pos_v = cols[pos_col].strip()
                    locus_key_chr = f"chr{chr_v}" if not chr_v.startswith("chr") else chr_v
                    locus_key = f"{locus_key_chr}:{pos_v}"
                else:
                    locus_key = None

                # Primary gene symbol
                primary = re.split(r"[,;]", gene_str)[0].strip().upper()
                if not primary or primary in (".", "-") or not primary.isalpha():
                    continue

                # Filtering: rare locus OR coding impact
                is_rare = (
                    locus_key is not None and
                    (locus_key in rare_loci or
                     locus_key in af_lookup)
                )
                is_coding = func_str in (
                    "exonic", "splicing", "exonic;splicing",
                    "ncRNA_exonic", "ncRNA_splicing",
                )

                if is_rare or is_coding:
                    genes[primary] = True

    except Exception as e:
        log.error("Failed to parse %s: %s", multianno_path, e)

    return sorted(genes.keys())


def run_single_vcf(args: Tuple[str, Path, float, int, Path, bool]) -> Dict:
    """
    Full pipeline for one VCF:
      1. Decompress .vcf.gz  (Python gzip)
      2. Find rare loci (AF < threshold + non-ref GT) + AF lookup
      3. Write pre-filtered VCF (rare loci only)
      4. convert2annovar.pl   → avinput  (on filtered VCF)
      5. table_annovar.pl     → refGene  (on filtered avinput)
      6. Parse multianno.txt   → gene symbols
      7. Write <case_key>_genes.txt
    """
    case_key, vcf_path, af_threshold, min_gq, out_dir, force = args
    _ensure_gzip()

    out_file = out_dir / f"{case_key}_genes.txt"
    if out_file.exists() and not force:
        return {"case": case_key, "status": "skipped", "n_genes": 0}

    tmp_base = out_dir / case_key
    tmp_vcf = out_dir / f"{case_key}_tmp.vcf"          # decompressed full VCF
    filtered_vcf = out_dir / f"{case_key}_filtered.vcf"  # rare loci only
    log_file = out_dir / f"{case_key}_run.log"

    perl_bin = str(STRAWBERRY_PERL).replace("\\", "/")
    annovar_scripts = str(ANNOVAR_DIR).replace("\\", "/")
    db_dir = str(ANNOVAR_DB).replace("\\", "/")
    filtered_str = str(filtered_vcf).replace("\\", "/")
    tmp_str = str(tmp_base).replace("\\", "/")

    # Build env with PERL5LIB so annotate_variation.pl (called by
    # table_annovar.pl via system()) can find Strawberry Perl modules
    perl_lib = str(ANNOVAR_DIR.parent / "Strawberry" / "perl" / "lib")
    import os as _os
    _env = dict(_os.environ)
    _env["PERL5LIB"] = perl_lib + ";" + _env.get("PERL5LIB", "")

    with open(log_file, "w", encoding="utf-8") as lfh:
        ts = time.strftime("%Y-%m-%d %H:%M:%S")
        lfh.write(f"[{ts}] VCF: {vcf_path}\n")

        # ── Step 1: Decompress ───────────────────────────────────────────────
        lfh.write(f"[{ts}] Step 1: Decompressing {vcf_path.name} ...\n")
        try:
            opener = _gzip_mod.open if vcf_path.suffix == ".gz" else open
            with opener(vcf_path, "rt", encoding="utf-8", errors="replace") as fin:
                with open(tmp_vcf, "w", encoding="utf-8") as fout:
                    for line in fin:
                        fout.write(line)
            lfh.write(f"[{ts}] Decompress OK: {tmp_vcf.stat().st_size:,} bytes\n")
        except Exception as e:
            lfh.write(f"[{ts}] Decompress FAILED: {e}\n")
            return {"case": case_key, "status": "error", "n_genes": 0}

        # ── Step 2: Find rare loci ──────────────────────────────────────────
        lfh.write(f"[{ts}] Step 2: Finding rare loci (AF < {af_threshold}, GT != 0/0) ...\n")
        try:
            af_lookup, rare_loci = _parse_vcf_find_rare_loci(
                tmp_vcf, af_threshold, min_gq
            )
            lfh.write(
                f"[{ts}] Rare loci: {len(rare_loci):,}  |  "
                f"AF entries: {len(af_lookup):,}\n"
            )
            if not rare_loci:
                lfh.write(f"[{ts}] No rare loci found — skipping ANNOVAR\n")
                # Still write an empty gene list
                with open(out_file, "w") as gf:
                    pass
                all_ok = True
                n_genes = 0
            else:
                all_ok = None  # placeholder
        except Exception as e:
            lfh.write(f"[{ts}] VCF parse FAILED: {e}\n")
            return {"case": case_key, "status": "error", "n_genes": 0}

        # ── Step 3: Write pre-filtered VCF ─────────────────────────────────
        if rare_loci:
            lfh.write(f"[{ts}] Step 3: Writing filtered VCF (rare loci only) ...\n")
            try:
                n_filtered = _write_filtered_vcf(tmp_vcf, rare_loci, filtered_vcf)
                lfh.write(f"[{ts}] Filtered VCF: {n_filtered:,} variants\n")
            except Exception as e:
                lfh.write(f"[{ts}] Filtered VCF FAILED: {e}\n")
                return {"case": case_key, "status": "error", "n_genes": 0}

            # ── Step 4: convert2annovar ──────────────────────────────────────
            cmd_conv = [
                perl_bin,
                f"{annovar_scripts}/convert2annovar.pl",
                filtered_str,
                "-format", "vcf",
                "-outfile", f"{tmp_str}.avinput",
                "-withzyg",
            ]
            lfh.write(f"\n[{ts}] Step 4: convert2annovar\n")
            lfh.write(f"[{ts}] CMD: {' '.join(cmd_conv)}\n")
            all_ok = True
            try:
                result = subprocess.run(
                    cmd_conv, cwd=str(ANNOVAR_DIR),
                    capture_output=True, text=True, timeout=600,
                    env=_env,
                )
                lfh.write(f"[{ts}] returncode={result.returncode}\n")
                if result.stdout:
                    lfh.write(f"[stdout]\n{result.stdout[:500]}\n")
                if result.stderr:
                    lfh.write(f"[stderr]\n{result.stderr[:500]}\n")
                if result.returncode != 0:
                    all_ok = False
            except subprocess.TimeoutExpired:
                lfh.write(f"[{ts}] TIMEOUT in convert2annovar\n")
                all_ok = False
            except Exception as e:
                lfh.write(f"[{ts}] Exception: {e}\n")
                all_ok = False

            # ── Step 5: table_annovar ────────────────────────────────────────
            multianno_file: Optional[Path] = None
            if all_ok:
                cmd_table = [
                    perl_bin,
                    f"{annovar_scripts}/table_annovar.pl",
                    f"{tmp_str}.avinput",
                    db_dir,
                    "-buildver", "hg19",
                    "-out", tmp_str,
                    "-protocol", "refGene",
                    "-operation", "g",
                    "-nastring", ".",
                ]
                lfh.write(f"\n[{ts}] Step 5: table_annovar\n")
                lfh.write(f"[{ts}] CMD: {' '.join(cmd_table)}\n")
                try:
                    result = subprocess.run(
                        cmd_table, cwd=str(ANNOVAR_DIR),
                        capture_output=True, text=True, timeout=600,
                        env=_env,
                    )
                    lfh.write(f"[{ts}] returncode={result.returncode}\n")
                    if result.stdout:
                        lfh.write(f"[stdout]\n{result.stdout[:500]}\n")
                    if result.stderr:
                        lfh.write(f"[stderr]\n{result.stderr[:500]}\n")
                    if result.returncode != 0:
                        all_ok = False
                    else:
                        for ext in (".hg19_multianno.txt", ".hg19_multianno.vcf"):
                            p = out_dir / f"{case_key}{ext}"
                            if p.exists():
                                multianno_file = p
                                break
                except subprocess.TimeoutExpired:
                    lfh.write(f"[{ts}] TIMEOUT in table_annovar\n")
                    all_ok = False
                except Exception as e:
                    lfh.write(f"[{ts}] Exception: {e}\n")
                    all_ok = False

            # ── Step 6: Parse → gene list ──────────────────────────────────
            if all_ok:
                if multianno_file is None:
                    multianno_file = out_dir / f"{case_key}.hg19_multianno.txt"
                if multianno_file and multianno_file.exists():
                    gene_list = extract_genes_from_multianno(
                        multianno_file, rare_loci, af_lookup, af_threshold
                    )
                    n_genes = len(gene_list)
                    with open(out_file, "w", encoding="utf-8") as gf:
                        for g in gene_list:
                            print(g, file=gf)
                    lfh.write(f"[{ts}] Extracted {n_genes} genes -> {out_file}\n")
                else:
                    all_ok = False
                    n_genes = 0
                    lfh.write(f"[{ts}] Multianno not found\n")
            else:
                n_genes = 0
        else:
            n_genes = 0

        # ── Cleanup ─────────────────────────────────────────────────────────
        for ext in (".avinput", ".avinput.vcf", ".hg19_multianno.txt",
                     ".hg19_multianno.vcf", ".hg19_refGene.txt",
                     ".refGene.variant_function", ".refGene.exonic_variant_function",
                     ".refGene.log", ".refGene.invalid_input",
                     ".refGene.fa", ".refGene.exonic_variant_function.orig"):
            try:
                f = out_dir / f"{case_key}{ext}"
                if f.exists():
                    f.unlink()
            except Exception:
                pass
        for f_path in (tmp_vcf, filtered_vcf):
            try:
                if f_path.exists():
                    f_path.unlink()
            except Exception:
                pass

    status = "ok" if all_ok else "error"
    return {"case": case_key, "status": status, "n_genes": n_genes}


# ── Dependency check ─────────────────────────────────────────────────────────

def check_dependencies() -> bool:
    errors = []
    if not STRAWBERRY_PERL.exists():
        errors.append(f"Strawberry Perl not found: {STRAWBERRY_PERL}")
    if not ANNOVAR_DIR.exists():
        errors.append(f"ANNOVAR not found: {ANNOVAR_DIR}")
    refgene = ANNOVAR_DB / "hg19_refGene.txt"
    if not refgene.exists():
        errors.append(f"ANNOVAR refGene missing: {refgene}")
    if errors:
        for e in errors:
            log.error("DEPENDENCY ERROR: %s", e)
        return False
    return True


# ── Main ─────────────────────────────────────────────────────────────────────

def discover_vcfs(case_filter: Optional[List[str]]) -> List[Tuple[str, Path]]:
    """Find all raw (non-VEP) VCF files matching HPO cases."""
    found: Dict[str, Path] = {}
    for vcf_path in VCF_DIR.glob("patient_case*_*.vcf.gz"):
        if vcf_path.name.endswith(".vep.vcf.gz"):
            continue
        case_id, sample_id = parse_case_id_from_vcf(vcf_path)
        case_key = f"{case_id}_{sample_id}"
        if case_filter and case_id not in case_filter:
            continue
        if case_key not in found:
            found[case_key] = vcf_path
    return sorted(found.items(), key=lambda x: x[0])


def main():
    parser = argparse.ArgumentParser(
        description="Extract candidate gene lists from raw VCFs using ANNOVAR.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Pipeline (per VCF):
  1. Decompress .vcf.gz (Python gzip)
  2. Parse INFO → AF; find rare loci (AF < 0.01 + non-ref GT)
  3. Write pre-filtered VCF (rare loci only)
  4. convert2annovar.pl  → avinput
  5. table_annovar.pl    → refGene annotation
  6. Parse multianno.txt → extract gene symbols
  7. Write <caseID>_<SampleID>_genes.txt

Examples:
  python extract_candidates.py                          # all cases
  python extract_candidates.py --workers 4             # parallel
  python extract_candidates.py --dry-run                # preview
  python extract_candidates.py --cases case1 case100   # specific cases
  python extract_candidates.py --af-threshold 0.001     # stricter rare filter
""",
    )
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite existing output files")
    parser.add_argument("--cases", nargs="+",
                        help="Specific case IDs (e.g. case1 case100)")
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument(
        "--af-threshold", type=float, default=DEFAULT_AF_THRESHOLD,
        help=f"Rare variant AF threshold (default: {DEFAULT_AF_THRESHOLD})",
    )
    parser.add_argument(
        "--min-gq", type=int, default=DEFAULT_MIN_GQ,
        help=f"Minimum genotype quality (GQ) (default: {DEFAULT_MIN_GQ})",
    )
    args = parser.parse_args()

    out_dir = args.out_dir if args.out_dir else OUT_DIR
    af_threshold = args.af_threshold
    min_gq = args.min_gq

    log.info("VCF dir : %s", VCF_DIR)
    log.info("Out dir : %s", out_dir)
    log.info("AF threshold: %.4f  |  Min GQ: %d", af_threshold, min_gq)

    log.info("Checking dependencies ...")
    if not check_dependencies():
        sys.exit(1)
    log.info("All dependencies OK.")

    out_dir.mkdir(parents=True, exist_ok=True)

    vcf_items = discover_vcfs(args.cases)
    log.info("Found %d VCF files to process.", len(vcf_items))

    already_done = sum(
        1 for case_key, _ in vcf_items
        if (out_dir / f"{case_key}_genes.txt").exists()
    )
    if args.force:
        to_process = len(vcf_items)
        already_done = 0
        log.info("Force mode: reprocessing all %d cases", to_process)
    else:
        to_process = len(vcf_items) - already_done
        log.info("Already complete: %d | To process: %d",
                 already_done, to_process)

    if args.dry_run:
        for case_key, vcf_path in vcf_items:
            out = out_dir / f"{case_key}_genes.txt"
            status = "skip" if out.exists() else "process"
            log.info("  [%s] %s  (%s)", status, case_key, vcf_path.name)
        sys.exit(0)

    results = []
    stats = {"ok": 0, "skipped": 0, "error": 0}

    if args.workers <= 1:
        for case_key, vcf_path in vcf_items:
            if not args.force and (out_dir / f"{case_key}_genes.txt").exists():
                stats["skipped"] += 1
                continue
            r = run_single_vcf((case_key, vcf_path, af_threshold, min_gq, out_dir, args.force))
            results.append(r)
            stats[r["status"]] += 1
            sym = "OK" if r["status"] == "ok" else r["status"]
            log.info("  %-30s %-8s  n_genes=%s", case_key, sym, r["n_genes"])
    else:
        work = [
            (case_key, vcf_path, af_threshold, min_gq, out_dir, args.force)
            for case_key, vcf_path in vcf_items
            if args.force or not (out_dir / f"{case_key}_genes.txt").exists()
        ]
        log.info("Starting parallel processing with %d workers ...", args.workers)
        with multiprocessing.Pool(processes=args.workers) as pool:
            for i, r in enumerate(pool.imap_unordered(run_single_vcf, work)):
                results.append(r)
                stats[r["status"]] += 1
                done = i + 1
                log.info("[%3d/%3d] ok=%-3d err=%-3d %-30s  n_genes=%s",
                         done, len(work), stats["ok"], stats["error"],
                         r["case"], r["n_genes"])

    log.info("========================================")
    log.info(" Done.")
    log.info("  total   = %d", len(vcf_items))
    log.info("  ok      = %d", stats["ok"])
    log.info("  skipped = %d", stats["skipped"])
    log.info("  errors  = %d", stats["error"])
    log.info(" Output   : %s", out_dir)
    log.info("========================================")


if __name__ == "__main__":
    main()
