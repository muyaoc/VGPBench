#!/usr/bin/env python3
"""
Unified Evaluation Script for All Gene-Prioritization Tools
============================================================
Reads the canonical master_table.csv (provided alongside this script) and
computes the full set of evaluation metrics reported in the manuscript:

  - N_submitted, N_processed, N_hit for each tool and Top-k threshold
  - Conditional Top-k accuracy   (N_hit / N_processed)
  - End-to-end Top-k accuracy    (N_hit / N_submitted, failures = miss)
  - Input compatibility           (N_processed / N_submitted)
  - Wilson 95% confidence intervals for every proportion
  - Failure-mode breakdown

The script then verifies every computed value against the reference numbers
embedded from Supplementary Tables S5 and S6, and reports PASS/FAIL.

master_table.csv is the single source of truth: one row per case (1,015
HPO-valid cases) with per-tool submitted/processed/failure-category/rank/
Top-k flags for all 8 tool configurations. It was built by parsing the raw
per-tool result files; this script turns it into the summary statistics and
confidence intervals reported in the manuscript.

Usage
-----
  python evaluate_all_tools.py                              # uses ./master_table.csv
  python evaluate_all_tools.py --master-table PATH.csv      # custom path
  python evaluate_all_tools.py --output-dir PATH            # where to write CSVs
  python evaluate_all_tools.py --no-verify                  # skip S5/S6 verification

Requirements: pandas, numpy.  (statsmodels optional for Wilson CI; a pure-Python
fallback is included so the script runs without it.)

Date: 2026-08-13
"""

import argparse
import math
import os
import sys
from collections import OrderedDict

import numpy as np
import pandas as pd

# ============================================================
# Defaults
# ============================================================
DEFAULT_MASTER_TABLE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'master_table.csv',
)
DEFAULT_OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output')

K_VALUES = [1, 5, 10, 20]
Z_95 = 1.959963984540054   # z for 95% CI

# ============================================================
# Tool configuration (matches build_master_table.py / manuscript)
# ============================================================
# Order preserved as in manuscript Table 1: HPO+VCF tools first, then HPO-only.
TOOLS = OrderedDict([
    ('AI-MARRVEL',      'HPO+VCF'),
    ('Exomiser',        'HPO+VCF'),
    ('ClinPrior-VCF',   'HPO+VCF'),
    ('ClinPrior-HPO',   'HPO-only'),
    ('Phenolyzer-VCF',  'HPO-only'),   # VCF defines candidate set; ranking is HPO-based
    ('Phenolyzer-Pheno','HPO-only'),
    ('Phen2Gene',       'HPO-only'),
    ('GADO',            'HPO-only'),
])

# Display labels used in the main-paper Table 1
DISPLAY_LABEL = {
    'AI-MARRVEL': 'AI-MARRVEL',
    'Exomiser': 'Exomiser',
    'ClinPrior-VCF': 'ClinPrior-VCF',
    'ClinPrior-HPO': 'ClinPrior-HPO',
    'Phenolyzer-VCF': 'Phenolyzer-VCFGenes',
    'Phenolyzer-Pheno': 'Phenolyzer-AllGenes',
    'Phen2Gene': 'Phen2Gene',
    'GADO': 'GADO',
}

# ============================================================
# Reference values from Supplementary Tables S5 & S6 (main.tex)
# Used for verification.  Format: tool -> {k: (n_sub, n_proc, n_hit, cond_pct, e2e_pct)}
# ============================================================
S5_REFERENCE = {
    'AI-MARRVEL': {
        1:  (777, 773, 260, 33.64, 33.46),
        5:  (777, 773, 460, 59.51, 59.20),
        10: (777, 773, 550, 71.15, 70.79),
        20: (777, 773, 618, 79.95, 79.54),
    },
    'Exomiser': {
        1:  (777, 777, 225, 28.96, 28.96),
        5:  (777, 777, 298, 38.35, 38.35),
        10: (777, 777, 330, 42.47, 42.47),
        20: (777, 777, 361, 46.46, 46.46),
    },
    'ClinPrior-VCF': {
        1:  (777, 265,  59, 22.26,  7.59),
        5:  (777, 265, 108, 40.75, 13.90),
        10: (777, 265, 136, 51.32, 17.50),
        20: (777, 265, 163, 61.51, 20.98),
    },
    'ClinPrior-HPO': {
        1:  (1015, 483, 20,  4.14, 1.97),
        5:  (1015, 483, 48,  9.94, 4.73),
        10: (1015, 483, 63, 13.04, 6.21),
        20: (1015, 483, 85, 17.60, 8.37),
    },
    'Phenolyzer-VCF': {
        1:  (777, 683, 17,  2.49, 2.19),
        5:  (777, 683, 40,  5.86, 5.15),
        10: (777, 683, 56,  8.20, 7.21),
        20: (777, 683, 78, 11.42, 10.04),
    },
    'Phenolyzer-Pheno': {
        1:  (1015, 548, 18,  3.28, 1.77),
        5:  (1015, 548, 52,  9.49, 5.12),
        10: (1015, 548, 69, 12.59, 6.80),
        20: (1015, 548, 89, 16.24, 8.77),
    },
    'Phen2Gene': {
        1:  (1015, 311,  6,  1.93, 0.59),
        5:  (1015, 311, 18,  5.79, 1.77),
        10: (1015, 311, 30,  9.65, 2.96),
        20: (1015, 311, 37, 11.90, 3.65),
    },
    'GADO': {
        1:  (1015, 909,  3,  0.33, 0.30),
        5:  (1015, 909,  6,  0.66, 0.59),
        10: (1015, 909, 11,  1.21, 1.08),
        20: (1015, 909, 14,  1.54, 1.38),
    },
}

# S6 compatibility reference: tool -> (n_sub, n_proc, n_failed, compat_pct, failure_pct)
S6_REFERENCE = {
    'AI-MARRVEL':       (777, 773,   4, 99.49,  0.51),
    'Exomiser':         (777, 777,   0, 100.00, 0.00),
    'ClinPrior-VCF':    (777, 265, 512, 34.11, 65.89),
    'ClinPrior-HPO':    (1015, 483, 532, 47.59, 52.41),
    'Phenolyzer-VCF':   (777, 683,  94, 87.90, 12.10),
    'Phenolyzer-Pheno': (1015, 548, 467, 53.99, 46.01),
    'Phen2Gene':        (1015, 311, 704, 30.64, 69.36),
    'GADO':             (1015, 909, 106, 89.56, 10.44),
}


# ============================================================
# Wilson score 95% confidence interval (pure Python, no statsmodels)
# ============================================================
def wilson_ci(k, n, z=Z_95):
    """Return (low, high) Wilson score interval for k successes in n trials."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1.0 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    spread = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, center - spread), min(1.0, center + spread))


def fmt_ci(low, high):
    """Format a (0-1) interval as 'LL.LL-HH.HH' percent."""
    return f"{low*100:.2f}-{high*100:.2f}"


# ============================================================
# Core metric computation
# ============================================================
def compute_metrics(df):
    """Compute per-tool, per-k metrics from master_table.

    Returns a list of dict rows (one per tool x k) and a compatibility summary.
    """
    rows = []
    compat_rows = []

    for tool, modality in TOOLS.items():
        sub_mask = df[f'{tool}_submitted'] == 1
        sub_df = df.loc[sub_mask]
        proc_df = df.loc[sub_mask & (df[f'{tool}_processed'] == 1)]

        n_sub = int(sub_mask.sum())
        n_proc = len(proc_df)
        n_failed = n_sub - n_proc

        # Compatibility
        compat_pct = 100.0 * n_proc / n_sub if n_sub else 0.0
        failure_pct = 100.0 * n_failed / n_sub if n_sub else 0.0
        c_low, c_high = wilson_ci(n_proc, n_sub)
        compat_rows.append({
            'Tool': DISPLAY_LABEL[tool],
            'Tool_key': tool,
            'Modality': modality,
            'N_submitted': n_sub,
            'N_processed': n_proc,
            'N_failed': n_failed,
            'Compatibility_pct': round(compat_pct, 2),
            'Compatibility_CI': fmt_ci(c_low, c_high),
            'Failure_pct': round(failure_pct, 2),
        })

        # Top-k metrics
        for k in K_VALUES:
            # Conditional: hits among processed
            n_hit_cond = int(proc_df[f'{tool}_top{k}'].sum())
            cond_pct = 100.0 * n_hit_cond / n_proc if n_proc else 0.0
            cond_low, cond_high = wilson_ci(n_hit_cond, n_proc)

            # End-to-end: hits among submitted (failures count as miss)
            n_hit_e2e = int(sub_df[f'{tool}_top{k}'].sum())
            e2e_pct = 100.0 * n_hit_e2e / n_sub if n_sub else 0.0
            e2e_low, e2e_high = wilson_ci(n_hit_e2e, n_sub)

            rows.append({
                'Tool': DISPLAY_LABEL[tool],
                'Tool_key': tool,
                'Modality': modality,
                'k': k,
                'N_submitted': n_sub,
                'N_processed': n_proc,
                'N_hit': n_hit_cond,
                'Conditional_Topk_pct': round(cond_pct, 2),
                'Conditional_CI': fmt_ci(cond_low, cond_high),
                'End_to_end_Topk_pct': round(e2e_pct, 2),
                'End_to_end_CI': fmt_ci(e2e_low, e2e_high),
            })

    return rows, compat_rows


# ============================================================
# Failure-mode breakdown
# ============================================================
def compute_failure_modes(df):
    """Return per-tool failure-category counts."""
    out = []
    for tool in TOOLS:
        sub_mask = df[f'{tool}_submitted'] == 1
        failed = df.loc[sub_mask & (df[f'{tool}_processed'] == 0)]
        if len(failed) == 0:
            out.append({'Tool': DISPLAY_LABEL[tool], 'Tool_key': tool,
                        'N_failed': 0, 'Failure_modes': '(none)'})
            continue
        cats = failed[f'{tool}_failure_category'].value_counts()
        cat_str = ', '.join(f'{cat}={cnt}' for cat, cnt in cats.items())
        out.append({'Tool': DISPLAY_LABEL[tool], 'Tool_key': tool,
                    'N_failed': int(len(failed)), 'Failure_modes': cat_str})
    return out


# ============================================================
# Verification against S5 / S6
# ============================================================
def approx_eq(a, b, tol=0.05):
    return abs(a - b) <= tol


def verify(rows, compat_rows):
    """Verify computed metrics against S5/S6 reference.  Returns (n_pass, n_fail, messages)."""
    msgs = []
    n_pass, n_fail = 0, 0

    # --- S5: Top-k accuracy ---
    row_lookup = {(r['Tool_key'], r['k']): r for r in rows}
    for tool, k_dict in S5_REFERENCE.items():
        for k, (ref_sub, ref_proc, ref_hit, ref_cond, ref_e2e) in k_dict.items():
            r = row_lookup.get((tool, k))
            if r is None:
                continue
            checks = [
                ('N_submitted', r['N_submitted'], ref_sub),
                ('N_processed', r['N_processed'], ref_proc),
                ('N_hit',       r['N_hit'],       ref_hit),
                ('Cond_pct',    r['Conditional_Topk_pct'], ref_cond),
                ('E2E_pct',     r['End_to_end_Topk_pct'],  ref_e2e),
            ]
            for label, got, ref in checks:
                ok = approx_eq(float(got), float(ref), tol=0.5)
                if ok:
                    n_pass += 1
                else:
                    n_fail += 1
                    msgs.append(f"  [S5] {tool} k={k} {label}: got {got}, ref {ref}  ✗")

    # --- S6: Compatibility ---
    compat_lookup = {r['Tool_key']: r for r in compat_rows}
    for tool, (ref_sub, ref_proc, ref_failed, ref_compat, ref_fail) in S6_REFERENCE.items():
        r = compat_lookup.get(tool)
        if r is None:
            continue
        checks = [
            ('N_submitted', r['N_submitted'], ref_sub),
            ('N_processed', r['N_processed'], ref_proc),
            ('N_failed',    r['N_failed'],    ref_failed),
            ('Compat_pct',  r['Compatibility_pct'], ref_compat),
            ('Failure_pct', r['Failure_pct'],       ref_fail),
        ]
        for label, got, ref in checks:
            ok = approx_eq(float(got), float(ref), tol=0.5)
            if ok:
                n_pass += 1
            else:
                n_fail += 1
                msgs.append(f"  [S6] {tool} {label}: got {got}, ref {ref}  ✗")

    return n_pass, n_fail, msgs


# ============================================================
# Pretty printing
# ============================================================
def print_summary_table(rows, compat_rows, failure_rows):
    """Print a manuscript-style summary table (Table 1 equivalent)."""
    print("\n" + "=" * 100)
    print("TABLE 1 EQUIVALENT — Tool Performance and Compatibility")
    print("=" * 100)
    hdr = (f"{'Modality':<10} {'Tool':<22} {'N_sub':>6} {'N_proc':>7} "
           f"{'Compat%':>8} {'CondT1%':>8} {'CondT5%':>8} {'CondT10%':>8} "
           f"{'CondT20%':>9} {'E2E_T1%':>8}")
    print(hdr)
    print("-" * 100)

    compat_by_tool = {r['Tool_key']: r for r in compat_rows}
    topk_by_tool = {}
    for r in rows:
        topk_by_tool.setdefault(r['Tool_key'], {})[r['k']] = r

    for tool, modality in TOOLS.items():
        c = compat_by_tool[tool]
        t = topk_by_tool[tool]
        print(f"{modality:<10} {c['Tool']:<22} {c['N_submitted']:>6} {c['N_processed']:>7} "
              f"{c['Compatibility_pct']:>8.2f} "
              f"{t[1]['Conditional_Topk_pct']:>8.2f} {t[5]['Conditional_Topk_pct']:>8.2f} "
              f"{t[10]['Conditional_Topk_pct']:>8.2f} {t[20]['Conditional_Topk_pct']:>9.2f} "
              f"{t[1]['End_to_end_Topk_pct']:>8.2f}")

    # Failure modes
    print("\n" + "=" * 100)
    print("FAILURE-MODE BREAKDOWN")
    print("=" * 100)
    for fr in failure_rows:
        print(f"  {fr['Tool']:<22} N_failed={fr['N_failed']:>4}  {fr['Failure_modes']}")


# ============================================================
# Main
# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description='Unified evaluation of all gene-prioritization tools from master_table.csv.')
    parser.add_argument('--master-table', default=DEFAULT_MASTER_TABLE,
                        help='Path to master_table.csv')
    parser.add_argument('--output-dir', default=DEFAULT_OUTPUT_DIR,
                        help='Directory to write summary CSVs')
    parser.add_argument('--no-verify', action='store_true',
                        help='Skip verification against S5/S6 reference values')
    args = parser.parse_args()

    master_path = args.master_table
    if not os.path.exists(master_path):
        print(f"ERROR: master_table.csv not found at {master_path}", file=sys.stderr)
        print(f"       Specify --master-table PATH to its location.", file=sys.stderr)
        sys.exit(1)

    print("=" * 100)
    print("Unified Tool Evaluation")
    print("=" * 100)
    print(f"Master table: {master_path}")

    df = pd.read_csv(master_path)
    print(f"Loaded {len(df)} cases, {len(df.columns)} columns")

    # Validate expected columns exist
    missing = [f'{t}_submitted' for t in TOOLS if f'{t}_submitted' not in df.columns]
    if missing:
        print(f"ERROR: master_table missing columns: {missing}", file=sys.stderr)
        sys.exit(1)

    # Compute metrics
    rows, compat_rows = compute_metrics(df)
    failure_rows = compute_failure_modes(df)

    # Print summary
    print_summary_table(rows, compat_rows, failure_rows)

    # Write CSVs
    os.makedirs(args.output_dir, exist_ok=True)
    s5_path = os.path.join(args.output_dir, 'evaluation_topk_with_ci.csv')
    s6_path = os.path.join(args.output_dir, 'evaluation_compatibility_with_ci.csv')
    fail_path = os.path.join(args.output_dir, 'evaluation_failure_modes.csv')

    pd.DataFrame(rows).drop(columns=['Tool_key']).to_csv(s5_path, index=False)
    pd.DataFrame(compat_rows).drop(columns=['Tool_key']).to_csv(s6_path, index=False)
    pd.DataFrame(failure_rows).drop(columns=['Tool_key']).to_csv(fail_path, index=False)

    print(f"\nCSV outputs written to {args.output_dir}/")
    print(f"  - evaluation_topk_with_ci.csv        (S5 equivalent)")
    print(f"  - evaluation_compatibility_with_ci.csv (S6 equivalent)")
    print(f"  - evaluation_failure_modes.csv        (S13 equivalent)")

    # Verification
    if not args.no_verify:
        print("\n" + "=" * 100)
        print("VERIFICATION against Supplementary Tables S5 & S6")
        print("=" * 100)
        n_pass, n_fail, msgs = verify(rows, compat_rows)
        total = n_pass + n_fail
        print(f"  Checks passed: {n_pass}/{total}")
        if n_fail == 0:
            print("\n  ✓ ALL METRICS MATCH S5/S6 REFERENCE VALUES")
        else:
            print(f"\n  ✗ {n_fail} MISMATCH(es) detected:")
            for m in msgs:
                print(m)
            print("\n  NOTE: Minor discrepancies may reflect paper-internal inconsistencies")
            print("        between main_paper.tex Table 1 and main.tex S5 (e.g. Exomiser Top-1")
            print("        28.83% vs 28.96%, ClinPrior-HPO 4.35% vs 4.14%). The master_table.csv")
            print("        is the canonical source of truth and matches S5 hit counts.")

    print("\nDone.")


if __name__ == '__main__':
    main()
