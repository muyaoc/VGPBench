# Composite Patient VCF Files

This directory will contain the composite virtual patient VCF files upon publication.

## Contents (to be released)

- Composite VCF files for virtual patients (1,041 cases): each file embeds a pathogenic variant within the assigned 1000 Genomes Phase 3 individual background, constructed using the pipeline in `scripts/virtual_patients/generate_patient_vcfs_bcftools_dedup.py`
- Tabix index files for each composite VCF
- VEP-annotated VCF files, required for ClinPrior-VCF mode

## Data will be made publicly available upon publication

The complete VCF dataset will be deposited in a permanent repository (e.g., Zenodo) and linked here upon manuscript publication. The actual file names and directory layout will be disclosed at that time.

## File format

- All VCFs are on GRCh37/hg19 coordinates
- Files are bgzip-compressed and tabix-indexed
