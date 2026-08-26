#!/usr/bin/env python3
"""
batch_phenolyzer.py

Batch wrapper for Phenolyzer — prioritizes candidate genes per case using
HPO phenotype terms and (optionally) a candidate gene list from VCF.

Supports three gene-source modes:
  vcf_genes   Extract protein-affecting gene symbols from VEP-annotated VCF
              and feed them as -gene candidates to Phenolyzer.
              This is the recommended mode for clinical WES/WGS analysis.
  all_genes   Use the full human gene universe as the candidate pool.
              Phenolyzer ranks all genes purely by HPO phenotype match.
              Useful as a baseline or when no VCF is available.
  phenotype   No gene candidate restriction — Phenolyzer generates a
              phenotype-only ranked list (shorthand for the original
              phenolyzer web-server behaviour).

Usage:
    python batch_phenolyzer.py --mode vcf_genes [--workers 4] [--dry-run]
    python batch_phenolyzer.py --mode all_genes [--workers 4]
    python batch_phenolyzer.py --mode phenotype [--workers 4]

Output:
    phenolyzer_results/
        case{N}_{sample}/
            case{N}_{sample}.seed_gene_list   ← ranked gene list (Phenolyzer output)
            case{N}_{sample}.final_gene_list  ← ranked + predicted genes
            case{N}_{sample}.merge_gene_scores
            case{N}_{sample}.annotated_gene_scores
        summary/
            phenolyzer_summary.tsv
            aggregated_genes.tsv   ← top-N genes across all cases

Requirements:
    - Perl with Bioperl (set PERL_BIN if not on PATH)
    - Parallel::ForkManager (installed via cpan)
    - Python 3.8+ (uses built-in gzip module — no external gzip needed)
"""

import os
import sys
import re
import csv
import time
import logging
import tempfile
import argparse
import subprocess
import multiprocessing
from pathlib import Path
from collections import defaultdict
from typing import Optional, List, Tuple, Dict, Any

# ── Logging ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)-7s %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("batch_phenolyzer")


# ── Paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent.resolve()
# Perl binary: set PERL_BIN if perl is not on PATH (e.g. on Windows:
#   set PERL_BIN=D:\Strawberry\perl\bin\perl.exe)
STRAWBERRY_PERL = Path(os.environ.get("PERL_BIN", "perl"))
PHENOLYZER_DIR = Path(os.environ.get(
    "PHENOLYZER_DIR", SCRIPT_DIR / "phenolyzer"
))
HPO_DIR = SCRIPT_DIR / "phenotype" / "GRCh37"
VCF_DIR = SCRIPT_DIR / "genotype" / "GRCh37"
OUT_DIR = SCRIPT_DIR / "phenolyzer_results"
HUMAN_GENELIST = PHENOLYZER_DIR / "lib" / "compiled_database" / "DB_HUMAN_GENE_ID"
CANDIDATE_DIR = SCRIPT_DIR / "phenolyzer_candidates"


# ── VEP CSQ parsing ──────────────────────────────────────────────────────────
# VEP VCF INFO/CSQ field format (each ALT allele has one CSQ entry, '|'-delimited):
#   0-Allele  1-Consequence  2-IMPACT  3-SYMBOL  4-Gene  5-Strand  6-VARIANT_CLASS
#   7-SYMBOL_SOURCE  8-HGNC_ID  ... (varies by VEP version)
#
# We extract SYMBOL (col 3) for variants with HIGH or MODERATE impact
# that overlap coding sequences.

IMPACTFUL_CONSEQUENCES = frozenset({
    # HIGH impact
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "protein_altering_variant",
    # MODERATE impact
    "missense_variant",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant_splicing",
    "splice_region_variant",
})

IMPACTFUL_SUBSTRINGS = ["missense", "nonsense", "frameshift", "splice",
                         "stop_", "start_", "inframe", "protein_altering"]


def _is_impactful(cs: str) -> bool:
    """Return True if the SO consequence term is protein-affecting."""
    if not cs:
        return False
    cs_lower = cs.lower()
    if cs_lower in IMPACTFUL_CONSEQUENCES:
        return True
    for kw in IMPACTFUL_SUBSTRINGS:
        if kw.lower() in cs_lower:
            return True
    return False


def extract_vcf_genes(vcf_path: str,
                      include_low: bool = False) -> List[str]:
    """
    Extract unique gene symbols from a VEP-annotated .vcf.gz (or .vcf).

    Returns:
        Sorted list of unique gene symbols (upper-case).
    """
    opener = _gzip_mod.open if vcf_path.endswith(".gz") else open
    genes: set = set()
    sym_idx = -1

    try:
        with opener(vcf_path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith("##INFO=<ID=CSQ"):
                    # Extract headers after "Format: " (handle quoted/descriptive headers)
                    idx = line.find("Format:")
                    if idx >= 0:
                        header_str = line[idx + 8:].strip()
                        # Remove leading " and trailing "
                        header_str = re.sub(r'^"', '', header_str)
                        header_str = re.sub(r'".*$', '', header_str)
                        headers = header_str.split("|")
                        for i, h in enumerate(headers):
                            if h == "SYMBOL":
                                sym_idx = i
                                break
                    break

            if sym_idx < 0:
                log.warning("SYMBOL column not found in CSQ header of %s", vcf_path)
                return []

            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) < 8:
                    continue
                info = fields[7]
                # CSQ may not be the first field in INFO
                csq_pos = info.find("CSQ=")
                if csq_pos < 0:
                    continue

                csq_str = info[csq_pos + 4:]
                for allele_entry in csq_str.split(","):
                    parts = allele_entry.split("|")
                    if len(parts) <= sym_idx:
                        continue
                    sym = parts[sym_idx].strip()
                    if not sym or not sym.isalpha():
                        continue
                    impact = parts[2].strip().upper() if len(parts) > 2 else ""
                    if impact in ("HIGH", "MODERATE") or _is_impactful(parts[1]):
                        genes.add(sym.upper())
                    elif include_low and impact == "LOW":
                        genes.add(sym.upper())
    except Exception as e:
        log.error("Failed to extract genes from VCF %s: %s", vcf_path, e)
        raise

    return sorted(genes)


# Lazy import so the script doesn't crash on a clean machine before gzip check
_gzip_mod = None  # type: ignore


def _ensure_gzip():
    global _gzip_mod
    if _gzip_mod is None:
        import gzip
        _gzip_mod = gzip


def extract_all_human_genes() -> List[str]:
    """Read all gene symbols from DB_HUMAN_GENE_ID."""
    _ensure_gzip()
    if not HUMAN_GENELIST.exists():
        log.error("Human gene list not found: %s", HUMAN_GENELIST)
        return []
    genes = []
    with open(HUMAN_GENELIST) as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)  # skip header
        for row in reader:
            if row and len(row) >= 2:
                genes.append(row[1].strip().upper())
    log.info("Loaded %d genes from human gene database.", len(genes))
    return sorted(genes)


def find_vcf(case_id: str, sample_id: str, vcf_dir: Path) -> Optional[Path]:
    """
    Find the VEP-annotated VCF for a given case.

    Priority: .vep.vcf.gz > .vcf.gz
    """
    candidates = [
        vcf_dir / f"patient_{case_id}_{sample_id}.vep.vcf.gz",
        vcf_dir / f"patient_{case_id}_{sample_id}.vcf.gz",
    ]
    for p in candidates:
        if p.exists():
            return p
    return None


def load_pregenerated_candidates(case_id: str, sample_id: str) -> Optional[List[str]]:
    """
    Read pre-generated candidate gene list from phenolyzer_candidates/.
    File format: <caseID>_<SampleID>_genes.txt  (one gene symbol per line)
    Returns None if the file doesn't exist.
    """
    gene_file = CANDIDATE_DIR / f"{case_id}_{sample_id}_genes.txt"
    if not gene_file.exists():
        return None
    genes = []
    with open(gene_file, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            gene = line.strip().upper()
            if gene:
                genes.append(gene)
    return genes if genes else None


def parse_hpo_file(hpo_file: Path) -> Tuple[List[str], int]:
    """Read HPO terms from a Phenolyzer-format HPO file (one HP: per line)."""
    terms = []
    with open(hpo_file) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("HP:"):
                terms.append(line)
    return terms, len(terms)


def parse_case_id(hpo_file: Path) -> Tuple[str, str]:
    """
    Parse case_id and sample_id from HPO filename like case1_HG00407_hpos.txt.

    Strips the '_hpos' suffix first, then splits on the last underscore
    to separate case ID from sample ID.
    """
    name = hpo_file.stem          # e.g. case1_HG00407_hpos
    name = re.sub(r'_hpos$', '', name)  # strip trailing _hpos
    parts = name.rsplit("_", 1)
    if len(parts) == 2:
        return parts[0], parts[1]
    return name, ""


# ── Phenolyzer runner ─────────────────────────────────────────────────────────

def build_phenolyzer_cmd(hpo_file: Path,
                          gene_list: Optional[List[str]],
                          out_prefix: Path,
                          nproc_phenolyzer: int,
                          use_precalc: bool,
                          use_prediction: bool = True) -> Tuple[List[str], Optional[str]]:
    """
    Build the full perl command for Phenolyzer.

    Returns:
        (cmd_list, gene_file_path) — caller must delete gene_file_path.
    """
    script = PHENOLYZER_DIR / "disease_annotation.pl"
    gene_file: Optional[str] = None

    # Use forward slashes throughout to avoid Perl regex issues on Windows
    perl_bin = str(STRAWBERRY_PERL).replace("\\", "/")
    script_path = str(script).replace("\\", "/")
    hpo_path = str(hpo_file).replace("\\", "/")

    cmd = [
        perl_bin,
        script_path,
        hpo_path,
        "-file",
        "-phenotype",
        "-logistic",
        "-out", str(out_prefix).replace("\\", "/"),
        "-addon", "DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE",
        "-addon_weight", "0.25",
        "-nproc", str(nproc_phenolyzer),
    ]

    if use_prediction:
        cmd.append("-prediction")

    if use_precalc and (PHENOLYZER_DIR / "lib" / "precalculated_expansions").exists():
        cmd.append("-use_precalc")

    mentha = PHENOLYZER_DIR / "lib" / "compiled_database" / "DB_MENTHA_GENE_GENE_INTERACTION"
    if mentha.exists():
        cmd += ["-addon_gg", "DB_MENTHA_GENE_GENE_INTERACTION",
                "-addon_gg_weight", "0.05"]

    # Write candidate gene list to a temp file if provided
    if gene_list:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False, encoding="utf-8"
        ) as gf:
            for g in gene_list:
                print(g, file=gf)
            gene_file = gf.name
        cmd += ["-gene", gene_file.replace("\\", "/")]

    return cmd, gene_file


def run_single_case(args: Tuple[str, Path, Optional[List[str]], str, bool, Path]) -> Dict[str, Any]:
    """
    Worker function for one case — runs Phenolyzer.

    args: (case_key, hpo_file, gene_list, gene_source, use_prediction, out_dir)
    Returns: dict with status info
    """
    case_key, hpo_file, gene_list, gene_source, use_prediction, out_dir = args
    _ensure_gzip()

    out_subdir = out_dir / case_key
    out_subdir.mkdir(parents=True, exist_ok=True)
    out_prefix = out_subdir / case_key
    log_file = out_subdir / "run.log"

    # Resume check — seed_gene_list exists regardless of prediction mode
    seed_list = str(out_subdir / f"{case_key}.seed_gene_list")
    if Path(seed_list).exists():
        return {
            "case": case_key,
            "status": "skipped",
            "reason": "already_complete",
            "gene_source": gene_source,
            "gene_count": len(gene_list) if gene_list else 0,
        }

    n_hpo, _ = parse_hpo_file(hpo_file)

    cmd, gene_file = build_phenolyzer_cmd(
        hpo_file=hpo_file,
        gene_list=gene_list,
        out_prefix=out_prefix,
        nproc_phenolyzer=1,          # internal Phenolyzer parallelism
        use_precalc=True,
        use_prediction=use_prediction,
    )

    # Convert all paths to forward slashes to avoid Perl regex issues on Windows
    cmd = [arg.replace("\\", "/") for arg in cmd]

    with open(log_file, "w", encoding="utf-8") as lfh:
        ts = time.strftime("%Y-%m-%d %H:%M:%S")
        lfh.write(f"[{ts}] {' '.join(cmd)}\n")
        lfh.write(f"[{ts}] gene_source={gene_source} n_genes={len(gene_list) if gene_list else 0}\n")
        lfh.write(f"[{ts}] n_hpo={n_hpo}\n")

        try:
            result = subprocess.run(
                cmd,
                cwd=None,  # run from current working directory so relative -out paths resolve correctly
                capture_output=True,
                text=True,
                timeout=600,  # 10 min per case
                input="A\n",  # "A" = overwrite all existing precalculated files
            )
            lfh.write(f"[{ts}] returncode={result.returncode}\n")
            if result.stderr:
                lfh.write(f"[stderr]\n{result.stderr[:2000]}\n")
            if result.stdout:
                lfh.write(f"[stdout]\n{result.stdout[:2000]}\n")

            ok = result.returncode == 0
            if ok:
                status = "ok"
                # seed_gene_list is always produced (with or without -prediction)
                seed_out = out_subdir / f"{case_key}.seed_gene_list"
                if not seed_out.exists():
                    ok = False
                    status = "error_no_output"
            else:
                status = "error_nonzero"

        except subprocess.TimeoutExpired:
            status = "error_timeout"
            lfh.write(f"[{ts}] TIMEOUT after 1800s\n")
            ok = False
        except Exception as e:
            status = "error_exception"
            lfh.write(f"[{ts}] Exception: {e}\n")
            ok = False

    # Clean up temp gene file
    if gene_file and os.path.exists(gene_file):
        try:
            os.unlink(gene_file)
        except Exception:
            pass

    return {
        "case": case_key,
        "status": status,
        "gene_source": gene_source,
        "gene_count": len(gene_list) if gene_list else 0,
        "n_hpo": n_hpo,
        "out_dir": str(out_subdir),
    }


# ── Summary aggregation ───────────────────────────────────────────────────────

def aggregate_top_genes(top_n: int = 100) -> Dict[str, float]:
    """
    Walk phenolyzer_results/ and collect top-N genes per case.
    Returns: dict of gene -> list of (case, rank, score)
    """
    gene_data: Dict[str, List[Tuple[str, int, float]]] = defaultdict(list)
    out_dir = OUT_DIR

    for case_dir in sorted(out_dir.iterdir()):
        if not case_dir.is_dir():
            continue
        seed_list = case_dir / f"{case_dir.name}.seed_gene_list"
        if not seed_list.exists():
            continue
        try:
            with open(seed_list, encoding="utf-8", errors="replace") as fh:
                reader = csv.reader(fh, delimiter="\t")
                header = next(reader, None)
                if header is None:
                    continue
                for row in reader:
                    if len(row) < 4:
                        continue
                    try:
                        rank = int(row[0])
                        gene = row[1].strip().upper()
                        score = float(row[3])
                    except (ValueError, IndexError):
                        continue
                    if rank > top_n:
                        break
                    gene_data[gene].append((case_dir.name, rank, score))
        except Exception as e:
            log.warning("Failed to parse %s: %s", seed_list, e)

    return dict(gene_data)


def write_summary(results: List[Dict[str, Any]], summary_tsv: Path):
    """Write per-case summary TSV."""
    with open(summary_tsv, "w", newline="", encoding="utf-8") as sfh:
        writer = csv.DictWriter(
            sfh,
            fieldnames=["case", "n_hpo", "n_genes", "gene_source", "status", "out_dir"],
            delimiter="\t",
        )
        writer.writeheader()
        for r in results:
            writer.writerow({
                "case": r["case"],
                "n_hpo": r.get("n_hpo", ""),
                "n_genes": r.get("gene_count", ""),
                "gene_source": r["gene_source"],
                "status": r["status"],
                "out_dir": r.get("out_dir", ""),
            })


def write_aggregated(gene_data: Dict[str, List[Tuple[str, int, float]]],
                     out_file: Path, top_cases: int = 5):
    """Write aggregated gene summary."""
    with open(out_file, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["Gene", "NumCases", "AvgRank", "TopRanks", "Cases"])
        for gene, appearances in sorted(
            gene_data.items(),
            key=lambda x: len(x[1]),  # sort by number of cases
            reverse=True,
        ):
            ranks = [r for _, r, _ in appearances]
            avg_rank = sum(ranks) / len(ranks)
            top_r = ",".join(str(r) for _, r, _ in sorted(appearances, key=lambda x: x[1])[:top_cases])
            cases = ",".join(c for c, _, _ in appearances[:top_cases])
            writer.writerow([gene, len(appearances), f"{avg_rank:.1f}", top_r, cases])


# ── Pre-run dependency checks ─────────────────────────────────────────────────

def check_dependencies() -> bool:
    """Verify all required tools are available."""
    errors = []

    if not STRAWBERRY_PERL.exists():
        errors.append(f"Strawberry Perl not found at: {STRAWBERRY_PERL}")
    else:
        # Test Parallel::ForkManager
        result = subprocess.run(
            [str(STRAWBERRY_PERL), "-e",
             "use Parallel::ForkManager; print 'OK'"],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            errors.append(
                "Parallel::ForkManager not installed in Perl. "
                "Run: perl -MCPAN -e 'CPAN::Shell->install(\"Parallel::ForkManager\")'"
            )

    script = PHENOLYZER_DIR / "disease_annotation.pl"
    if not script.exists():
        errors.append(f"Phenolyzer disease_annotation.pl not found at: {script}")

    if not HUMAN_GENELIST.exists():
        errors.append(f"Human gene list not found at: {HUMAN_GENELIST}")

    # Note: Python's gzip module (import gzip) is used — no external gzip binary needed.

    if errors:
        for e in errors:
            log.error("DEPENDENCY ERROR: %s", e)
        return False
    return True


# ── Case discovery ────────────────────────────────────────────────────────────

def discover_cases(hpo_dir: Path,
                   mode: str,
                   vcf_dir: Path,
                   all_genes: Optional[List[str]]) -> List[Tuple[str, Path, Optional[List[str]], str, bool]]:
    """
    Walk HPO files and build the full case list with gene candidates.

    Returns:
        List of (case_key, hpo_file, gene_list_or_None, gene_source, use_prediction)
    """
    hpo_files = sorted(hpo_dir.glob("*_hpos.txt"))
    if not hpo_files:
        log.error("No *_hpos.txt files found in %s", hpo_dir)
        sys.exit(1)

    log.info("Found %d HPO files in %s", len(hpo_files), hpo_dir)

    cases: List[Tuple[str, Path, Optional[List[str]], str, bool]] = []

    for i, hpo_file in enumerate(hpo_files):
        case_id, sample_id = parse_case_id(hpo_file)
        case_key = f"{case_id}_{sample_id}"

        gene_list: Optional[List[str]] = None
        gene_source = ""
        use_prediction = True  # default

        if mode == "vcf_genes":
            # Read from pre-generated candidate gene lists (by extract_candidates.py)
            gene_list = load_pregenerated_candidates(case_id, sample_id)
            if gene_list:
                gene_source = f"vcf_pregen:{CANDIDATE_DIR.name}/{case_id}_{sample_id}_genes.txt ({len(gene_list)} genes)"
            else:
                log.warning("No pre-generated candidate list for %s — falling back to human_all", case_key)
                if all_genes is None:
                    _ensure_gzip()
                    all_genes = extract_all_human_genes()
                gene_list = all_genes[:]
                gene_source = "human_all_fallback"

            use_prediction = True

        elif mode == "all_genes":
            # No gene candidate restriction — use the full Phenolyzer database.
            # -prediction triggers O(N^2) gene-gene loops that crash on Windows;
            # the DB already covers all associations so ranking is identical.
            gene_list = None
            gene_source = "phenolyzer_db_all"
            use_prediction = False

        elif mode == "phenotype":
            gene_list = None
            gene_source = "none"
            use_prediction = False

        cases.append((case_key, hpo_file, gene_list, gene_source, use_prediction))

        if (i + 1) % 100 == 0:
            log.info("  ... processed %d/%d cases", i + 1, len(hpo_files))

    return cases


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Batch Phenolyzer runner for HPO-based gene prioritization.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Modes:
  vcf_genes   Extract coding-gene candidates from VEP-annotated VCF and feed
              them as -gene seeds to Phenolyzer. (default)
  all_genes   Use the full human gene universe (~20 000 genes) as the candidate
              pool. Phenolyzer ranks all genes purely by HPO phenotype match.
              Baseline / fallback when no VCF is available.
  phenotype   No gene candidate restriction — phenotype-only ranking.
              Equivalent to the original Phenolyzer web-server behaviour.

Examples:
  # Recommended: VCF-derived gene candidates, 4 parallel workers
  python batch_phenolyzer.py --mode vcf_genes --workers 4

  # Baseline: rank all genes by phenotype only
  python batch_phenolyzer.py --mode all_genes --workers 4

  # Dry-run to preview what would be run
  python batch_phenolyzer.py --mode vcf_genes --dry-run
""",
    )
    parser.add_argument(
        "--mode", choices=["vcf_genes", "all_genes", "phenotype"],
        default="vcf_genes",
        help="Gene source mode (default: vcf_genes)",
    )
    parser.add_argument(
        "--workers", type=int, default=1,
        help="Number of parallel case workers (default: 1)",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print case list without executing Phenolyzer",
    )
    parser.add_argument(
        "--hpo-dir", default=None,
        help="HPO files directory (default: ./phenotype/GRCh37)",
    )
    parser.add_argument(
        "--vcf-dir", default=None,
        help="VCF files directory (default: ./genotype/GRCh37)",
    )
    parser.add_argument(
        "--out-dir", default=None,
        help="Output directory (default: ./phenolyzer_results)",
    )
    parser.add_argument(
        "--perl-bin", default=None,
        help="Path to Perl binary (default: perl on PATH; set PERL_BIN or this flag if different)",
    )
    parser.add_argument(
        "--include-low-impact", action="store_true",
        help="Also include LOW-impact (synonymous) variants as gene candidates "
             "(not recommended — adds noise)",
    )
    parser.add_argument(
        "--top-n", type=int, default=100,
        help="Top-N genes per case for aggregated output (default: 100)",
    )
    args = parser.parse_args()

    # ── Apply path overrides ───────────────────────────────────────────────────
    # These globals are read by discover_cases(), run_single_case(), etc.
    global HPO_DIR, VCF_DIR, OUT_DIR, STRAWBERRY_PERL
    if args.hpo_dir is not None:
        HPO_DIR = Path(args.hpo_dir)
    if args.vcf_dir is not None:
        VCF_DIR = Path(args.vcf_dir)
    if args.out_dir is not None:
        OUT_DIR = Path(args.out_dir)
    if args.perl_bin is not None:
        STRAWBERRY_PERL = Path(args.perl_bin)

    log.info("Mode: %s | Workers: %d | Dry-run: %s", args.mode, args.workers, args.dry_run)

    log.info("HPO dir : %s", HPO_DIR)
    log.info("VCF dir : %s", VCF_DIR)
    log.info("Out dir : %s", OUT_DIR)

    # ── Dependency check ─────────────────────────────────────────────────────
    log.info("Checking dependencies ...")
    if not check_dependencies():
        sys.exit(1)
    log.info("All dependencies OK.")

    # ── Patch PATH for subprocess workers ──────────────────────────────────────
    # On Windows, Perl needs unzip.exe (bundled with Git for Windows) for
    # precalculated expansions. Set GIT_BIN to Git's usr/bin, or leave it unset
    # on systems where unzip is already on PATH.
    git_bin = os.environ.get("GIT_BIN", "")
    if git_bin and os.path.isdir(git_bin):
        os.environ["PATH"] = git_bin + os.pathsep + os.environ.get("PATH", "")
        log.info("Patched PATH: %s prepended (for unzip).", git_bin)

    # ── Setup output ─────────────────────────────────────────────────────────
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_dir = OUT_DIR / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    summary_tsv = summary_dir / "phenolyzer_summary.tsv"
    aggregated_tsv = summary_dir / "aggregated_genes.tsv"

    # ── Discover cases ───────────────────────────────────────────────────────
    all_genes: Optional[List[str]] = None
    cases = discover_cases(
        hpo_dir=HPO_DIR,
        mode=args.mode,
        vcf_dir=VCF_DIR,
        all_genes=all_genes,
    )

    # ── Pre-flight: check how many would be skipped ───────────────────────────
    total = len(cases)
    def is_complete(case_key: str) -> bool:
        return (OUT_DIR / case_key / f"{case_key}.seed_gene_list").exists()

    already_done = sum(1 for case_key, *_ in cases if is_complete(case_key))
    log.info("Total cases: %d | Already complete: %d | To process: %d",
             total, already_done, total - already_done)

    if args.dry_run:
        log.info("[DRY RUN] Would process %d cases:", total - already_done)
        for case_key, hpo_file, gene_list, gene_source, _pred_flag in cases:
            if is_complete(case_key):
                continue
            log.info("  %s  source=%s  genes=%d", case_key, gene_source,
                     len(gene_list) if gene_list else 0)
        sys.exit(0)

    # ── Run cases ─────────────────────────────────────────────────────────────
    # Prediction is expensive (O(N^2) gene-gene network loops) — skip for
    # all_genes/phenotype modes since the DB already covers all associations
    use_prediction = (args.mode == "vcf_genes")

    # Build work items (only non-complete cases)
    work_items = [
        (case_key, hpo_file, gene_list, gene_source, use_prediction, OUT_DIR)
        for case_key, hpo_file, gene_list, gene_source, _pred_flag in cases
        if not is_complete(case_key)
    ]

    results: List[Dict[str, Any]] = []
    stats = defaultdict(int)

    if args.workers <= 1:
        # Serial processing
        for item in work_items:
            r = run_single_case(item)
            results.append(r)
            stats[r["status"]] += 1
            status_sym = "OK" if r["status"] == "ok" else r["status"]
            log.info("[%3d/%3d] %-30s %s", stats["ok"] + stats.get("error_nonzero", 0) +
                      stats.get("error_timeout", 0) + stats.get("skipped", 0) + 1,
                      len(work_items), r["case"], status_sym)
    else:
        # Parallel processing via multiprocessing.Pool
        log.info("Starting parallel processing with %d workers ...", args.workers)
        with multiprocessing.Pool(processes=args.workers) as pool:
            for i, r in enumerate(pool.imap_unordered(run_single_case, work_items)):
                results.append(r)
                stats[r["status"]] += 1
                done = i + 1
                ok = stats["ok"]
                err = sum(v for k, v in stats.items() if k.startswith("error"))
                log.info("[%3d/%3d] ok=%-3d err=%-3d %s", done, len(work_items),
                         ok, err, r["case"])

    # ── Write summaries ───────────────────────────────────────────────────────
    # Merge skipped (already-complete) cases into results for summary
    for case_key, hpo_file, gene_list, gene_source, _pred_flag in cases:
        if is_complete(case_key):
            results.append({
                "case": case_key,
                "status": "skipped",
                "gene_source": gene_source,
                "gene_count": len(gene_list) if gene_list else 0,
                "n_hpo": len(parse_hpo_file(hpo_file)[0]),
                "out_dir": str(OUT_DIR / case_key),
            })

    write_summary(results, summary_tsv)
    log.info("Summary written to: %s", summary_tsv)

    # ── Aggregated gene summary ──────────────────────────────────────────────
    log.info("Building aggregated gene summary ...")
    gene_data = aggregate_top_genes(top_n=args.top_n)
    write_aggregated(gene_data, aggregated_tsv)
    log.info("Aggregated genes written to: %s", aggregated_tsv)

    # ── Final report ─────────────────────────────────────────────────────────
    log.info("========================================")
    log.info(" Done.")
    log.info("  total  = %d", total)
    log.info("  ok     = %d", stats["ok"])
    log.info("  skipped= %d", stats["skipped"])
    for k, v in sorted(stats.items()):
        if k not in ("ok", "skipped"):
            log.info("  %-8s= %d", k, v)
    log.info(" Results : %s", OUT_DIR)
    log.info(" Summary : %s", summary_tsv)
    log.info("========================================")


if __name__ == "__main__":
    main()
