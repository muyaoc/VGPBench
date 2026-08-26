#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Batch convert HGVS notations from published_diagnosed_pathogenic_variants.merged.txt
to per-case VCF files using VEP Docker for reproducibility.

Compared to hgvs_to_vcf_ensembl.py (Ensembl REST API), this script provides:
  - Reproducibility: Docker image version is pinned -> VEP version + DB version are fixed
  - Version traceability: VCF header records full VEP/Ensembl version info

How it works:
  VEP's --cache mode (without --offline) can resolve HGVS input:
  - HGVS resolution: connects to Ensembl DB (version matched to Docker image)
  - Variant annotation: uses local cache (fast)

Requirements:
  - Docker (with ensemblorg/ensembl-vep image pulled)
  - VEP cache directory (set via --vep-cache-dir or GENETOOLS_BASE env var)

Usage:
  # Basic usage
  python hgvs_to_vcf_vep_docker.py \\
    --merged published_diagnosed_pathogenic_variants_combined.txt \\
    --outdir pathogenic_variants_docker

  # Test a single case
  python hgvs_to_vcf_vep_docker.py \\
    --merged published_diagnosed_pathogenic_variants_combined.txt \\
    --outdir pathogenic_variants_docker \\
    --case 1

  # Download GRCh38 cache (needed before processing GRCh38 cases)
  python hgvs_to_vcf_vep_docker.py --download-grch38-cache
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import OrderedDict, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


# ============================================================
# Default configuration
# ============================================================
DEFAULT_DOCKER_IMAGE = "ensemblorg/ensembl-vep:latest"
# VEP cache dir: derive from GENETOOLS_BASE env var, or pass --vep-cache-dir.
_env_base = os.environ.get("GENETOOLS_BASE")
DEFAULT_VEP_CACHE_DIR = os.path.join(_env_base, "vep_data") if _env_base else None
ASSEMBLY_MAP = {
    "GRCh37": "GRCh37",
    "GRCh38": "GRCh38",
    "hg19": "GRCh37",
    "hg38": "GRCh38",
}


# ============================================================
# Argument parsing
# ============================================================
def parse_args():
    p = argparse.ArgumentParser(
        description="Convert HGVS to per-case VCF using pinned VEP Docker image "
                    "for reproducibility. "
                    "Replaces hgvs_to_vcf_ensembl.py (REST API) with VEP Docker."
    )
    p.add_argument(
        "--merged",
        default="published_diagnosed_pathogenic_variants_combined.txt",
        help="Path to merged pathogenic table (TSV). Default: %(default)s",
    )
    p.add_argument(
        "--outdir",
        default="pathogenic_variants_docker",
        help="Base output dir for per-case VCFs. Default: %(default)s",
    )
    p.add_argument(
        "--docker-image",
        default=DEFAULT_DOCKER_IMAGE,
        help=f"VEP Docker image tag. Default: %(default)s",
    )
    p.add_argument(
        "--vep-cache-dir",
        default=DEFAULT_VEP_CACHE_DIR,
        help="Local VEP cache directory (mounted into Docker). Default: %(default)s",
    )
    p.add_argument(
        "--case",
        help="Process only this specific case ID (e.g. '1')",
    )
    p.add_argument(
        "--cases",
        nargs="+",
        help="Process specific case IDs (e.g. 1 2 3)",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output VCF files",
    )
    p.add_argument(
        "--dedup-info",
        action="store_true",
        help="When duplicate coordinates are found, append all HGVS notations "
             "to the INFO field. Without this flag, only the first HGVS is kept.",
    )
    p.add_argument(
        "--download-grch38-cache",
        action="store_true",
        help="Download GRCh38 VEP cache and exit. "
             "Required before processing GRCh38 cases for the first time.",
    )
    p.add_argument(
        "--keep-vep-output",
        action="store_true",
        help="Keep the raw VEP output VCF (with CSQ annotations) alongside "
             "the simplified per-case VCF.",
    )
    p.add_argument(
        "--timeout",
        type=int,
        default=300,
        help="Docker run timeout in seconds. Default: %(default)s",
    )
    return p.parse_args()


# ============================================================
# HGVS loading (same logic as hgvs_to_vcf_ensembl.py)
# ============================================================
def load_hgvs_by_case(merged_path: str) -> Dict[Tuple[str, str], List[str]]:
    """
    Read the merged table, group HGVS by (Case, Reference_genome).

    Column structure:
      #Case   Gene   cHGVS   pHGVS   Ref seq   Phenotype   Disease   Reference_genome
    """
    if not os.path.isfile(merged_path):
        sys.stderr.write(f"[ERROR] merged file not found: {merged_path}\n")
        sys.exit(1)

    by_case: Dict[Tuple[str, str], List[str]] = defaultdict(list)

    with open(merged_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            case = str(row["#Case"]).strip()
            genome = row["Reference_genome"].strip()
            refseq = row["Ref seq"].strip()
            chgvs_raw = (row.get("cHGVS") or "").strip()

            if not case or not genome or not refseq:
                continue
            if not chgvs_raw or chgvs_raw == "-":
                continue

            # One row may have multiple cHGVS, comma-separated
            for c in chgvs_raw.split(","):
                c = c.strip()
                if not c or c == "-":
                    continue
                hgvs = f"{refseq}:{c}"
                by_case[(case, genome)].append(hgvs)

    # Deduplicate by HGVS string
    for key in list(by_case.keys()):
        uniq = sorted(set(by_case[key]))
        by_case[key] = uniq

    return by_case


# ============================================================
# Docker / VEP utility functions
# ============================================================
def check_docker() -> str:
    """Check Docker availability and return version string."""
    try:
        result = subprocess.run(
            ["docker", "--version"],
            capture_output=True, text=True, check=True, timeout=10
        )
        return result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        return ""


def check_docker_image(image: str) -> bool:
    """Check if the VEP Docker image is available locally."""
    try:
        result = subprocess.run(
            ["docker", "images", "-q", image],
            capture_output=True, text=True, check=True, timeout=10
        )
        return bool(result.stdout.strip())
    except Exception:
        return False


def check_vep_cache(cache_dir: str, assembly: str) -> bool:
    """Check if VEP cache exists for the given assembly."""
    cache_path = Path(cache_dir) / "homo_sapiens"
    if not cache_path.exists():
        return False
    for d in cache_path.iterdir():
        if d.is_dir() and assembly in d.name:
            return True
    return False


def get_vep_cache_version(cache_dir: str, assembly: str) -> Optional[str]:
    """Get the VEP cache version for the given assembly."""
    cache_path = Path(cache_dir) / "homo_sapiens"
    if not cache_path.exists():
        return None
    for d in cache_path.iterdir():
        if d.is_dir() and assembly in d.name:
            # e.g., "115_GRCh37" -> "115"
            parts = d.name.split("_")
            if parts:
                return parts[0]
    return None


def download_grch38_cache(cache_dir: str, docker_image: str):
    """Download GRCh38 VEP cache using VEP Docker's INSTALL.pl."""
    print(f"[INFO] Downloading GRCh38 VEP cache...")
    print(f"[INFO] Cache directory: {cache_dir}")
    print(f"[INFO] Docker image: {docker_image}")
    print(f"[INFO] This may take 10-30 minutes depending on network speed.")
    print()

    cmd = [
        "docker", "run", "--rm",
        "-v", f"{cache_dir}:/opt/vep/.vep",
        docker_image,
        "perl", "INSTALL.pl",
        "-a", "c",       # cache only
        "-s", "homo_sapiens",
        "-y",            # auto-yes
        "-g", "GRCh38",
    ]

    print(f"[CMD] {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False, text=True)
    if result.returncode != 0:
        sys.stderr.write("[ERROR] Cache download failed.\n")
        sys.exit(1)

    print(f"\n[OK] GRCh38 cache downloaded to {cache_dir}")


# ============================================================
# VEP Docker execution
# ============================================================
def run_vep_docker(
    hgvs_list: List[str],
    assembly: str,
    docker_image: str,
    vep_cache_dir: str,
    output_dir: str,
    timeout: int = 300,
) -> Tuple[str, str, int]:
    """
    Run VEP Docker to resolve HGVS notations.

    Args:
        hgvs_list:     List of HGVS notations
        assembly:      GRCh37 or GRCh38
        docker_image:  Docker image tag
        vep_cache_dir: Local VEP cache directory
        output_dir:    Local directory for Docker I/O
        timeout:       Docker run timeout in seconds

    Returns:
        (stdout, stderr, return_code)
    """
    # Create input file
    input_file = os.path.join(output_dir, "vep_input.txt")
    with open(input_file, "w") as f:
        for hgvs in hgvs_list:
            f.write(hgvs + "\n")

    cmd = [
        "docker", "run", "--rm",
        "-v", f"{vep_cache_dir}:/opt/vep/.vep:ro",
        "-v", f"{input_file}:/opt/vep/input.txt:ro",
        "-v", f"{output_dir}:/opt/vep/output",
        docker_image,
        "vep",
        "--cache",
        "--dir_cache", "/opt/vep/.vep",
        "--species", "homo_sapiens",
        "--assembly", assembly,
        "--format", "hgvs",
        "-i", "/opt/vep/input.txt",
        "-o", "/opt/vep/output/vep_output.vcf",
        "--vcf",
        "--no_stats",
        "--force_overwrite",
    ]

    result = subprocess.run(
        cmd, capture_output=True, text=True, timeout=timeout
    )

    return result.stdout, result.stderr, result.returncode


def parse_vep_vcf_output(vcf_path: str) -> Tuple[List[dict], str]:
    """
    Parse VEP VCF output to extract variant records and VEP version info.

    Returns:
        (records, vep_version_string)
        records: list of dicts with keys: chrom, pos, ref, alt, hgvs
    """
    records = []
    vep_version = ""

    if not os.path.exists(vcf_path):
        return records, vep_version

    with open(vcf_path, "r") as f:
        for line in f:
            if line.startswith("##VEP="):
                # Extract VEP version string
                vep_version = line.strip()
            elif line.startswith("#"):
                continue
            else:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 5:
                    continue

                chrom = cols[0]
                pos = cols[1]
                hgvs_id = cols[2]  # VEP puts the HGVS notation in the ID column
                ref = cols[3]
                alt = cols[4]

                # Handle multiple ALT alleles (comma-separated)
                alt_alleles = alt.split(",")
                for a in alt_alleles:
                    rec = {
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": a.strip(),
                        "hgvs": hgvs_id.replace(" ", "_"),
                    }
                    records.append(rec)

    return records, vep_version


def resolve_hgvs_robust(
    hgvs_list: List[str],
    assembly: str,
    docker_image: str,
    vep_cache_dir: str,
    timeout: int = 300,
) -> Tuple[List[dict], List[Tuple[str, str]], str]:
    """
    Resolve HGVS notations using VEP Docker, with fallback for failures.

    Strategy:
      1. Try all HGVS at once (fastest, one Docker run)
      2. If VEP crashes, try each HGVS individually to isolate failures
      3. Return (successful_records, failed_items, vep_version)

    Args:
        hgvs_list:     List of HGVS notations
        assembly:      GRCh37 or GRCh38
        docker_image:  Docker image tag
        vep_cache_dir: Local VEP cache directory
        timeout:       Docker run timeout per run

    Returns:
        (records, failures, vep_version)
        records: list of dicts with keys: chrom, pos, ref, alt, hgvs
        failures: list of (hgvs, reason) tuples
    """
    with tempfile.TemporaryDirectory(prefix="vep_docker_") as tmpdir:
        # Make output dir writable by Docker
        output_subdir = os.path.join(tmpdir, "output")
        os.makedirs(output_subdir, exist_ok=True)
        os.chmod(output_subdir, 0o777)

        # Step 1: Try all HGVS at once
        stdout, stderr, rc = run_vep_docker(
            hgvs_list, assembly, docker_image, vep_cache_dir,
            output_subdir, timeout=timeout
        )

        vcf_path = os.path.join(output_subdir, "vep_output.vcf")

        if rc == 0 and os.path.exists(vcf_path) and os.path.getsize(vcf_path) > 0:
            records, vep_version = parse_vep_vcf_output(vcf_path)

            # Check if all HGVS were resolved
            resolved_hgvs = {r["hgvs"] for r in records}
            failures = []
            for h in hgvs_list:
                if h.replace(" ", "_") not in resolved_hgvs:
                    failures.append((h, "resolved but no VCF record output"))

            return records, failures, vep_version

        # Step 2: VEP crashed -- try each HGVS individually
        print(f"    [FALLBACK] VEP batch failed (exit code {rc}), trying individually...")

        all_records = []
        all_failures = []
        vep_version = ""

        for hgvs in hgvs_list:
            single_dir = tempfile.mkdtemp(prefix="vep_single_")
            single_output = os.path.join(single_dir, "output")
            os.makedirs(single_output, exist_ok=True)
            os.chmod(single_output, 0o777)

            try:
                sout, serr, src = run_vep_docker(
                    [hgvs], assembly, docker_image, vep_cache_dir,
                    single_output, timeout=timeout
                )

                single_vcf = os.path.join(single_output, "vep_output.vcf")

                if src == 0 and os.path.exists(single_vcf) and os.path.getsize(single_vcf) > 0:
                    recs, ver = parse_vep_vcf_output(single_vcf)
                    if recs:
                        all_records.extend(recs)
                        if ver:
                            vep_version = ver
                    else:
                        all_failures.append((hgvs, "no VCF record in output"))
                else:
                    # Parse error from stderr
                    err_msg = "unknown error"
                    for line in serr.strip().split("\n"):
                        if "MSG:" in line:
                            err_msg = line.strip()
                            break
                    all_failures.append((hgvs, err_msg[:200]))
            except subprocess.TimeoutExpired:
                all_failures.append((hgvs, "timeout"))
            except Exception as e:
                all_failures.append((hgvs, str(e)[:200]))
            finally:
                try:
                    shutil.rmtree(single_dir)
                except Exception:
                    pass

        return all_records, all_failures, vep_version


# ============================================================
# Coordinate deduplication
# ============================================================
def dedup_by_coordinate(
    records: List[dict], keep_all_hgvs: bool = False
) -> List[dict]:
    """
    Deduplicate by (CHROM, POS, REF, ALT).
    """
    seen = OrderedDict()
    dup_count = 0

    for rec in records:
        coord = (rec["chrom"], rec["pos"], rec["ref"], rec["alt"])
        if coord not in seen:
            seen[coord] = rec
        else:
            dup_count += 1
            if keep_all_hgvs:
                existing = seen[coord]
                if "extra_hgvs" not in existing:
                    existing["extra_hgvs"] = []
                existing["extra_hgvs"].append(rec["hgvs"])

    if dup_count > 0:
        print(f"    [DEDUP] Removed {dup_count} duplicate coordinate(s)")

    return list(seen.values())


# ============================================================
# VCF writing
# ============================================================
def write_case_vcf(
    case: str,
    genome: str,
    records: List[dict],
    outdir: str,
    vep_version: str = "",
    keep_vep_output: bool = False,
    vep_raw_vcf: Optional[str] = None,
):
    """
    Write a simplified per-case VCF file (compatible with hgvs_to_vcf_ensembl.py output format).
    """
    subdir = os.path.join(outdir, genome)
    os.makedirs(subdir, exist_ok=True)
    out_path = os.path.join(subdir, f"case{case}.vcf")

    with open(out_path, "w", encoding="utf-8") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##source=VEP_Docker\n")
        out.write(f"##reference_genome={genome}\n")
        if vep_version:
            out.write(f"{vep_version}\n")
        out.write('##INFO=<ID=HGVS,Number=1,Type=String,Description="Original HGVS notation">\n')
        out.write('##INFO=<ID=HGVS_ALT,Number=.,Type=String,Description="Additional HGVS notations at same coordinate (deduped)">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for rec in records:
            chrom = rec["chrom"]
            pos = rec["pos"]
            ref = rec["ref"]
            alt = rec["alt"]
            hgvs = rec["hgvs"]
            info_fields = [f"HGVS={hgvs}"]

            if "extra_hgvs" in rec and rec["extra_hgvs"]:
                alt_hgvs = ",".join(rec["extra_hgvs"])
                info_fields.append(f"HGVS_ALT={alt_hgvs}")

            info = ";".join(info_fields)
            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{info}\n")

    print(f"[OK] Case {case} ({genome}), variants: {len(records)}, VCF: {out_path}")

    # Optionally copy raw VEP output
    if keep_vep_output and vep_raw_vcf and os.path.exists(vep_raw_vcf):
        raw_out = os.path.join(subdir, f"case{case}_vep_raw.vcf")
        shutil.copy2(vep_raw_vcf, raw_out)
        print(f"    Raw VEP output saved to: {raw_out}")


# ============================================================
# Main function
# ============================================================
def main():
    args = parse_args()

    # VEP cache dir must be provided (via --vep-cache-dir or GENETOOLS_BASE env var)
    if not args.vep_cache_dir:
        sys.stderr.write(
            "[ERROR] VEP cache directory not specified.\n"
            "        Set GENETOOLS_BASE (e.g. export GENETOOLS_BASE=/path/to/genetools)\n"
            "        or pass --vep-cache-dir /path/to/vep_cache\n"
        )
        sys.exit(1)

    # Handle --download-grch38-cache
    if args.download_grch38_cache:
        download_grch38_cache(args.vep_cache_dir, args.docker_image)
        return

    # Check Docker
    docker_version = check_docker()
    if not docker_version:
        sys.stderr.write("[ERROR] Docker is not installed or not in PATH\n")
        sys.exit(1)
    print(f"[INFO] Docker: {docker_version}")

    # Check Docker image
    if not check_docker_image(args.docker_image):
        sys.stderr.write(
            f"[ERROR] Docker image '{args.docker_image}' not found locally.\n"
            f"  Pull it first: docker pull {args.docker_image}\n"
        )
        sys.exit(1)
    print(f"[INFO] VEP Docker image: {args.docker_image}")

    # Check VEP cache
    for assembly in ["GRCh37", "GRCh38"]:
        has_cache = check_vep_cache(args.vep_cache_dir, assembly)
        cache_ver = get_vep_cache_version(args.vep_cache_dir, assembly)
        status = f"v{cache_ver}" if has_cache else "NOT FOUND"
        print(f"[INFO] VEP cache ({assembly}): {status}")
    print()

    # Load HGVS by case
    merged_path = os.path.abspath(args.merged)
    outdir = os.path.abspath(args.outdir)
    print(f"[INFO] Merged table : {merged_path}")
    print(f"[INFO] Output dir   : {outdir}")

    by_case = load_hgvs_by_case(merged_path)
    print(f"[INFO] Total (Case,Genome) groups: {len(by_case)}")

    # Filter by --case
    if args.case:
        target_case = str(args.case).strip()
        by_case = {k: v for k, v in by_case.items() if k[0] == target_case}
        print(f"[INFO] Filtering for Case: {target_case}")

    # Filter by --cases
    if args.cases:
        case_set = set(args.cases)
        by_case = {k: v for k, v in by_case.items() if k[0] in case_set}
        print(f"[INFO] Filtering for Cases: {args.cases}")

    if not by_case:
        print("[WARN] No cases to process after filtering.")
        return

    # Process each case
    total_success = 0
    total_failed = 0
    total_skipped = 0
    total_dedup = 0
    total_hgvs_failures = 0

    for (case, genome), hgvs_list in sorted(
        by_case.items(), key=lambda x: (int(x[0][0]) if x[0][0].isdigit() else 999, x[0][1])
    ):
        # Map genome name
        assembly = ASSEMBLY_MAP.get(genome, genome)

        # Check if output already exists
        out_vcf = os.path.join(outdir, assembly, f"case{case}.vcf")
        if os.path.exists(out_vcf) and not args.overwrite:
            print(f"[SKIP] Case {case} ({assembly}): output already exists")
            total_skipped += 1
            continue

        # Check cache availability
        if not check_vep_cache(args.vep_cache_dir, assembly):
            if assembly == "GRCh38":
                print(f"[WARN] Case {case}: No VEP cache for {assembly}, skipping. "
                      f"Use --download-grch38-cache to download.")
            else:
                print(f"[WARN] Case {case}: No VEP cache for {assembly}, skipping.")
            total_failed += 1
            continue

        print(f"\n{'='*50}")
        print(f"Case {case} / {assembly}")
        print(f"  HGVS count: {len(hgvs_list)}")
        start_time = time.time()

        # Run VEP Docker
        records, failures, vep_version = resolve_hgvs_robust(
            hgvs_list, assembly, args.docker_image, args.vep_cache_dir,
            timeout=args.timeout
        )

        elapsed = time.time() - start_time

        # Report failures
        if failures:
            print(f"  [FAIL] {len(failures)} HGVS could not be resolved:")
            for h, reason in failures:
                print(f"    - {h}: {reason[:100]}")
            total_hgvs_failures += len(failures)

        if not records:
            print(f"  [ERROR] No valid records for Case {case} ({assembly})")
            total_failed += 1
            continue

        # Deduplicate by coordinate
        before_dedup = len(records)
        records = dedup_by_coordinate(records, keep_all_hgvs=args.dedup_info)
        after_dedup = len(records)
        dup_count = before_dedup - after_dedup
        total_dedup += dup_count
        if dup_count > 0:
            print(f"  [DEDUP] {before_dedup} -> {after_dedup} records "
                  f"({dup_count} duplicate(s) removed by coordinate)")

        # Write per-case VCF
        write_case_vcf(
            case, assembly, records, outdir,
            vep_version=vep_version,
            keep_vep_output=args.keep_vep_output,
        )

        print(f"  Time: {elapsed:.1f}s")
        total_success += 1

    # Summary
    print(f"\n{'='*60}")
    print(f"Processing complete!")
    print(f"  Successful:         {total_success}")
    print(f"  Failed:             {total_failed}")
    print(f"  Skipped:            {total_skipped}")
    print(f"  HGVS failures:      {total_hgvs_failures}")
    print(f"  Duplicates removed: {total_dedup}")
    print(f"  Docker image:       {args.docker_image}")
    print(f"  Output directory:   {outdir}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
