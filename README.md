# VGPBench: Virtual paired Genotype and Phenotype Benchmark

This repository contains the scripts, reference data, and documentation for constructing and evaluating **VGPBench**, a standardized paired genotype-and-phenotype benchmark dataset for phenotype-driven gene prioritization.

**Associated manuscript:** A link to the associated manuscript will be provided here upon formal publication.

---

## Overview

This benchmark integrates molecularly diagnosed monogenic disease cases from published literature with population variant backgrounds from the 1000 Genomes Project to create virtual patients with explicit causal-gene labels. Unlike benchmarks based solely on disease-level annotations, this dataset preserves patient-level phenotypic sparsity and heterogeneity while providing standardized inputs in both HPO (Human Phenotype Ontology) and VCF (Variant Call Format) formats.

**Key statistics:**
- 1,289 curated literature-derived cases
- 1,015 virtual patients with valid HPO annotations
- 777 dual-modality cases (both HPO + VCF)
- 238 HPO-only cases
- 264 VCF-only cases
- 761 unique genes, 408 unique diseases, 20 clinical specialties

---

## Repository Structure

```
├── README.md                          # This file
├── LICENSE                            # License information
├── data/
│   ├── reference/                     # Reference data
│   │   ├── README.md                                                # ⏳ Master variant table placeholder (1,289 cases, released upon publication)
│   │   └── phenotype_hpo_cache.json                               # JAX HPO API query cache (141 entries, for reproducibility)
│   ├── genotype/                     # ⏳ Composite VCF files (will be made publicly available upon publication)
│   │   └── GRCh37/                    #   1,041 composite VCFs (bgzip + tabix indexed)
│   └── phenotype/                    # ⏳ HPO annotation files (will be made publicly available upon publication)
│       ├── GRCh37/                    #   779 HPO files for paired HPO+VCF cases
│       └── no_chgvs_cases/           #   236 HPO files for HPO-only cases
├── scripts/
│   ├── virtual_patients/              # Virtual patient construction pipeline
│   │   ├── hgvs_to_vcf_ensembl_dedup.py     # HGVS-to-VCF conversion via Ensembl VEP REST API
│   │   ├── hgvs_to_vcf_vep_docker.py         # HGVS-to-VCF conversion via local VEP Docker
│   │   ├── hgvs_to_vcf_38_19.py              # GRCh38-to-GRCh37 liftOver for VCF coordinates
│   │   ├── generate_patient_vcfs_bcftools_dedup.py  # Composite VCF assembly (bcftools merge)
│   │   └── map_phenotypes_to_hpo_improved.py # Phenotype-to-HPO standardization and validation
│   ├── tools/                         # Tool execution scripts
│   │   ├── AI-MARRVEL/
│   │   │   ├── run_aim.sh                     # AI-MARRVEL batch execution script
│   │   │   ├── nextflow.config                # Nextflow configuration (Docker profile)
│   │   │   ├── main.nf                        # Nextflow workflow definition
│   │   │   └── conf/                          # Nextflow config modules (base, modules)
│   │   ├── exomiser/
│   │   │   └── auto_exomiser_all_update.py    # Exomiser batch execution and configuration
│   │   ├── ClinPrior/
│   │   │   ├── run_clinprior_batch_new.R      # ClinPrior HPO+VCF batch (variant prioritization)
│   │   │   └── run_clinprior_batch_hpo.R      # ClinPrior HPO-only batch (gene prioritization → clinprior_results_hpo)
│   │   ├── phenolyzer/
│   │   │   ├── batch_phenolyzer.py            # Phenolyzer batch execution wrapper
│   │   │   └── extract_candidates.py          # VCF-to-candidate-gene list extraction
│   │   ├── Phen2Gene/
│   │   │   └── phen2gene-run.py               # Phen2Gene batch execution wrapper
│   │   └── GADO/
│   │       └── sort_GADO.py                   # GADO API query and result sorting
│   ├── evaluation/                    # Evaluation and metric calculation
│   │   ├── evaluate_all_tools.py              # ★ Unified evaluation: all 8 tools, Top-k + CIs + S5/S6 verification
│   │   └── master_table.csv                   # ★ Canonical per-case results (1,015 cases × 8 tools) — input to evaluate_all_tools.py
│   └── vep_annotation/                # VEP annotation preprocessing
│       ├── run_vep_single.sh                  # Single-case VEP annotation
│       └── run_vep_batch.sh                   # Batch VEP annotation
```

---

## Data Availability

### Currently Available

The following reference data are included in this repository:

| File | Description |
|------|-------------|
| `data/reference/phenotype_hpo_cache.json` | JAX HPO API query cache (141 phenotype→HPO mappings), generated during phenotype standardization for reproducibility. |
| `scripts/evaluation/master_table.csv` | Canonical per-case results table (1,015 cases × 8 tool configurations): submitted/processed/failure/rank/Top-k flags. Single source of truth for all statistical analyses; input to `evaluate_all_tools.py`. |

### Will Be Made Publicly Available Upon Publication

The following data will be released through a permanent repository (e.g., Zenodo) upon publication:

| Data | Description | Size (approx.) |
|------|-------------|----------------|
| `data/reference/published_diagnosed_pathogenic_variants_combined.txt` | Master variant table with 1,289 cases: causal gene, cHGVS, phenotype, disease, reference genome | ~200 KB |
| `data/genotype/GRCh37/` | 1,041 composite VCF files (bgzip + tabix indexed) | ~125 GB |
| `data/genotype/GRCh37/` (VEP-annotated) | VEP-annotated VCF files for ClinPrior | ~416 GB |
| `data/phenotype/GRCh37/` | 779 HPO annotation files (one HPO ID per line) | ~3 MB |
| `data/phenotype/no_chgvs_cases/` | 236 HPO-only annotation files | ~1 MB |

> **Note:** Composite VCF files embed pathogenic variants within 1000 Genomes Phase 3 population backgrounds. All patient-level phenotype and variant information is derived from published, publicly accessible literature and OMIM. No restricted or identifiable patient data are included.

---

## Pipeline: Virtual Patient Construction

The construction pipeline follows parallel genotype and phenotype branches:

### Genotype Branch

1. **HGVS-to-VCF conversion** (`hgvs_to_vcf_ensembl_dedup.py` or `hgvs_to_vcf_vep_docker.py`):
   - Convert cHGVS variants to genomic coordinates using the Ensembl VEP REST API or local VEP Docker
   - GRCh38 variants are lifted over to GRCh37 (`hgvs_to_vcf_38_19.py`)
   - Deduplication and reference-allele validation

2. **Composite VCF assembly** (`generate_patient_vcfs_bcftools_dedup.py`):
   - Normalize pathogenic variant VCFs
   - Remove overlapping background sites at pathogenic variant positions
   - Merge with assigned 1000 Genomes individual background VCF using `bcftools concat`
   - Zygosity-aware GT assignment based on inheritance mode

### Phenotype Branch

3. **Phenotype-to-HPO mapping** (`map_phenotypes_to_hpo_improved.py`):
   - Three-stage mapping: (1) extract embedded HP:XXXXXXX IDs, (2) match against local dictionary (79 common terms), (3) query JAX HPO API for remaining terms
   - HPO ID validation against a fixed release to remove obsolete/malformed identifiers
   - API query results cached for reproducibility

### Cohort Assignment

4. Cases with both valid HPO and VCF → **HPO+VCF cohort** (777 cases)
5. HPO-valid records without VCF → **HPO-only subset** (238 cases)
6. VCF-only records → **VCF-only subset** (264 cases)

---

## Evaluated Tools

| Tool | Version | Input modality | Execution |
|------|---------|----------------|-----------|
| AI-MARRVEL | 1.1.3 | HPO+VCF | Docker / Nextflow DSL2 |
| Exomiser | 14.0.0 | HPO+VCF | Local JVM (Java ≥17) |
| ClinPrior | 2.0 | HPO-only / HPO+VCF | Docker |
| Phenolyzer | — | HPO-only | Local Perl + Python |
| Phen2Gene | — | HPO-only | Local Python |
| GADO | API v1 | HPO-only | Remote REST API |

### Tool-specific notes

- **AI-MARRVEL**: Accepts raw VCF + HPO directly. Requires Docker with `zhandongliulab/aim-lite:1.2` image.
- **Exomiser**: Requires Exomiser data release 2406 for hg19. Optional CADD v1.4 and REMM data.
- **ClinPrior-VCF**: Requires VEP-annotated VCF input (see `scripts/vep_annotation/`). Minimum two HPO terms required.
- **ClinPrior-HPO**: HPO-only gene prioritization via `proteinScore()` + `MatrixPropagation()` (`run_clinprior_batch_hpo.R`). Minimum two HPO terms required; single-term cases fail at `proteinScore()` and are skipped (explaining the low input-compatibility rate).
- **Phenolyzer**: Two modes — AllGenes (all human genes as candidates) and VCFGenes (genes present in VCF only).
- **Phen2Gene**: May produce empty outputs under sparse HPO annotation (single HPO term).
- **GADO**: Remote API endpoint; HTTP 500 errors observed for some submissions.

---

## Evaluation pipeline

All manuscript statistics are reproducible from a single script plus the provided `master_table.csv`:

**Unified evaluation** (`master_table.csv` → summary metrics + verification):
```bash
python scripts/evaluation/evaluate_all_tools.py
```
Reads `scripts/evaluation/master_table.csv` (auto-detected next to the script) and computes conditional/end-to-end Top-k accuracy, Wilson 95% CIs, input compatibility, and failure-mode breakdown for all 8 tool configurations. Outputs CSVs matching Supplementary Tables S5, S6, and S13. Automatically verifies every computed value (N_submitted, N_processed, N_hit, accuracy %, compatibility %) against the S5/S6 reference values embedded from the manuscript — **200/200 checks pass**.

`master_table.csv` is the canonical per-case results table: one row per case (1,015 HPO-valid cases) with per-tool submitted/processed/failure-category/rank/Top-k flags for all 8 tool configurations. It was built by parsing the raw per-tool result files (Exomiser `.genes.tsv`, AI-MARRVEL `_gene_candidates.tsv`, ClinPrior CSVs, Phenolyzer gene lists, GADO result CSVs, Phen2Gene output files) and is provided here as the single source of truth for all statistical analyses.

---

## Requirements

### Software dependencies

- Python ≥3.8 (pandas, numpy, requests)
- R ≥3.6 (for ClinPrior)
- bcftools ≥1.20
- tabix
- Java ≥17 (for Exomiser)
- Docker (for AI-MARRVEL, ClinPrior)
- Nextflow ≥23.10.1 (for AI-MARRVEL)
- VEP (Ensembl Variant Effect Predictor)
- Perl with Bioperl (for Phenolyzer)

### Reference data

- 1000 Genomes Project Phase 3 VCFs (GRCh37/hg19)
- HPO ontology (version used in benchmark: see `map_phenotypes_to_hpo_improved.py`)
- Exomiser data release 2406
- AI-MARRVEL data dependencies (`aim-data-dependencies-2.4-public`)
- ClinPrior data v2.0-assemblyGRCh37

### Path configuration

The tool-execution and VEP-annotation scripts read input/output locations from the `GENETOOLS_BASE` environment variable rather than hard-coded paths. Set it to your local project root (where `genotype/`, `phenotype/`, and tool data directories live) before running:

```bash
export GENETOOLS_BASE=/path/to/genetools
```

Scripts that honor `GENETOOLS_BASE`:
- `scripts/tools/AI-MARRVEL/run_aim.sh`
- `scripts/tools/ClinPrior/run_clinprior_batch_new.R`
- `scripts/tools/ClinPrior/run_clinprior_batch_hpo.R`
- `scripts/tools/exomiser/auto_exomiser_all_update.py`
- `scripts/vep_annotation/run_vep_single.sh`
- `scripts/virtual_patients/hgvs_to_vcf_vep_docker.py` (or pass `--vep-cache-dir`)

For Phenolyzer, set `PERL_BIN` if `perl` is not on `PATH` (e.g. `export PERL_BIN=/path/to/perl`); on Windows you may also set `GIT_BIN` to Git's `usr/bin` so `unzip` is available to Perl.

---

## Reproducibility Notes

- All HPO IDs were validated against a fixed HPO release to ensure ontology updates do not alter benchmark inputs across runs.
- JAX HPO API query results are cached in `phenotype_hpo_cache.json` for reproducibility.
- The 1000 Genomes sample assignment (case → individual ID) will be provided as a mapping table alongside the VCF/HPO data upon publication.
- Composite VCFs use GRCh37/hg19 coordinates throughout; GRCh38 variants are lifted over before VCF assembly.

---

## Citation

If you use VGPBench or its data, please cite the associated manuscript. The citation and a link to the published paper will be added here upon formal publication.

> [Citation / manuscript link — to be added upon publication]

---

## License

This repository is released under the [MIT License](LICENSE). All data are derived from publicly accessible sources (1000 Genomes Project, OMIM, published literature) in compliance with applicable data usage agreements.

---

## Contact

For questions about the dataset or scripts, please open an issue on this repository or contact [corresponding author email TBD].
