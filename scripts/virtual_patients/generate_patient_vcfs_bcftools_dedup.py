#!/usr/bin/env python3
"""
Batch generate virtual patient VCF files using bcftools for proper merging.

This version uses bcftools to ensure proper VCF format compatibility with AI-MARRVEL.

Compared to generate_patient_vcfs_bcftools.py, this version adds:
  1. Deduplication of variant lines by (CHROM, POS, REF, ALT) in normalize_pathogenic_vcf
  2. **Pre-merge overlap removal**: before merging, exclude from the background VCF
     any sites that share the same (CHROM, POS) with the pathogenic variant VCF,
     so that the pathogenic variant always takes priority at overlapping positions.
  3. Summary reporting of duplicates found and background sites excluded.
  4. **Zygosity-aware GT assignment**: reads inheritance mode and zygosity from
     a .with_zygosity.txt data file to set the correct GT per variant:
       - AD  → 0/1 (heterozygous)
       - AR single-variant → 1/1 (homozygous)
       - AR compound het → 0/1 (each variant heterozygous)
       - X-linked male → 1 (hemizygous)
       - Mitochondrial → 1 (homoplasmy)
"""

import os
import sys
import subprocess
import argparse
import csv
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional


def read_sample_mapping(mapping_file: str) -> Dict[str, str]:
    """Read the case-to-sample mapping from used_samples_combine.txt"""
    mapping = {}
    with open(mapping_file, 'r') as f:
        header = f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                case_num = parts[0]
                sample_id = parts[1]
                mapping[case_num] = sample_id
    return mapping


def check_bcftools():
    """Check if bcftools is available"""
    try:
        result = subprocess.run(['bcftools', '--version'],
                              capture_output=True, text=True, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def read_zygosity_data(data_file: str) -> Dict[str, Dict]:
    """
    Read the zygosity-annotated variants data file.

    Returns:
        Dict mapping case_num (str) → {
            'inheritance_mode': str,
            'zygosity': str,
            'gt': str,
            'n_variants': int,
        }
    """
    data = {}
    if not os.path.exists(data_file):
        return data

    with open(data_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            case_num = row.get('#Case', '').strip()
            if not case_num:
                continue
            data[case_num] = {
                'inheritance_mode': row.get('Inheritance_mode', '').strip(),
                'zygosity': row.get('Zygosity', '').strip(),
                'gt': row.get('GT', '0/1').strip(),
                'n_variants': len([v for v in row.get('cHGVS', '').split(',') if v.strip()]) if row.get('cHGVS', '').strip() else 0,
            }
    return data


def normalize_pathogenic_vcf(input_vcf: str, output_vcf: str, sample_id: str,
                              gt_value: str = "0/1") -> Tuple[int, int]:
    """
    Normalize pathogenic variant VCF to match 1000 Genomes format
    - Add proper VCF headers including FORMAT definition
    - Add proper QUAL and FILTER fields
    - Ensure consistent genotype format
    - Sort variants by position
    - Deduplicate by (CHROM, POS, REF, ALT)
    - Set GT based on inheritance-aware zygosity (gt_value parameter)

    Args:
        input_vcf:  path to the raw pathogenic variant VCF
        output_vcf: path to write the normalized VCF
        sample_id:  sample identifier for the VCF column header
        gt_value:   genotype string to assign (e.g., "0/1", "1/1", "1")

    Returns:
        (total_variants, unique_variants) — counts before and after dedup
    """
    header_lines = []
    variant_lines = []
    has_format_header = False
    has_filter_header = False

    # Read and process the file
    with open(input_vcf, 'r') as f_in:
        for line in f_in:
            if line.startswith('##'):
                header_lines.append(line.rstrip('\n'))
                if '##FORMAT=<ID=GT' in line:
                    has_format_header = True
                if '##FILTER=<ID=PASS' in line:
                    has_filter_header = True
            elif line.startswith('#CHROM'):
                # Store column header for later
                col_header = line.rstrip('\n').split('\t')
            else:
                # Store variant lines for sorting
                cols = line.rstrip('\n').split('\t')
                if len(cols) >= 5:
                    variant_lines.append(cols)

    # Sort variants by chromosome and position
    def get_sort_key(cols):
        chrom = cols[0]
        try:
            pos = int(cols[1])
        except:
            pos = 0

        # Chromosome sorting: 1-22, X, Y, MT
        if chrom.isdigit():
            chrom_num = int(chrom)
        elif chrom == 'X':
            chrom_num = 23
        elif chrom == 'Y':
            chrom_num = 24
        elif chrom in ['MT', 'M']:
            chrom_num = 25
        else:
            chrom_num = 99

        return (chrom_num, pos)

    variant_lines.sort(key=get_sort_key)

    total_variants = len(variant_lines)

    # Deduplicate by (CHROM, POS, REF, ALT)
    seen_coords = set()
    deduped_lines = []
    for cols in variant_lines:
        coord_key = (cols[0], cols[1], cols[3], cols[4])  # CHROM, POS, REF, ALT
        if coord_key not in seen_coords:
            seen_coords.add(coord_key)
            deduped_lines.append(cols)

    dup_count = total_variants - len(deduped_lines)
    if dup_count > 0:
        print(f"    [DEDUP] Variant VCF: {total_variants} → {len(deduped_lines)} records "
              f"({dup_count} duplicate(s) removed by coordinate)")
    variant_lines = deduped_lines
    unique_variants = len(variant_lines)

    # Write normalized VCF
    with open(output_vcf, 'w') as f_out:
        # Write existing headers
        for header in header_lines:
            f_out.write(header + '\n')

        # Add missing FORMAT header if needed
        if not has_format_header:
            f_out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')

        # Add missing FILTER header if needed
        if not has_filter_header:
            f_out.write('##FILTER=<ID=PASS,Description="All filters passed">\n')

        # Write column header
        if len(col_header) < 10:
            # Add FORMAT and sample columns
            f_out.write('\t'.join(col_header) + '\tFORMAT\t' + sample_id + '\n')
        else:
            f_out.write('\t'.join(col_header) + '\n')

        # Write sorted & deduped variants
        for cols in variant_lines:
            # Ensure QUAL and FILTER are set
            if len(cols) > 5 and (cols[5] == '.' or cols[5] == ''):
                cols[5] = '100'
            if len(cols) > 6 and (cols[6] == '.' or cols[6] == ''):
                cols[6] = 'PASS'

            # Add FORMAT and genotype if missing
            if len(cols) < 9:
                cols.append('GT')
            if len(cols) < 10:
                cols.append(gt_value)

            f_out.write('\t'.join(cols) + '\n')

    return total_variants, unique_variants


def extract_pathogenic_positions(variant_vcf: str) -> Set[Tuple[str, int]]:
    """
    Read a pathogenic variant VCF and return the set of (chrom, pos) positions.

    These positions will be used to exclude overlapping sites from the background VCF.
    """
    positions = set()
    with open(variant_vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) >= 2:
                chrom = cols[0]
                try:
                    pos = int(cols[1])
                except ValueError:
                    continue
                positions.add((chrom, pos))
    return positions


def write_exclusion_bed(positions: Set[Tuple[str, int]], bed_path: str) -> int:
    """
    Write pathogenic variant positions to a BED file for bcftools region exclusion.

    BED format uses 0-based, half-open intervals:
      chrom  start  end
      1      12344  12345   (represents VCF position 12345)

    Returns:
        Number of positions written
    """
    # Sort for deterministic output
    sorted_pos = sorted(positions, key=lambda x: (
        int(x[0]) if x[0].isdigit() else (23 if x[0] == 'X' else 24 if x[0] == 'Y' else 25),
        x[1]
    ))

    with open(bed_path, 'w') as f:
        for chrom, pos in sorted_pos:
            # BED is 0-based, half-open: [pos-1, pos)
            f.write(f"{chrom}\t{pos - 1}\t{pos}\n")

    return len(sorted_pos)


def filter_background_exclude_positions(
    background_vcf: str,
    exclusion_bed: str,
    output_vcf: str,
    temp_dir: Path
) -> int:
    """
    Remove from the background VCF any sites whose (CHROM, POS) overlaps with
    the pathogenic variant positions, using bcftools view -T ^bed.gz.

    This ensures that at overlapping positions, the pathogenic variant is kept
    and the background record is removed.

    Args:
        background_vcf: path to the extracted single-sample background VCF
        exclusion_bed:  path to the BED file listing positions to exclude
        output_vcf:     path to write the filtered background VCF
        temp_dir:       temporary directory for intermediate files

    Returns:
        Number of background records removed at overlapping positions
    """
    # Count variants in original background
    def count_variants(vcf_path):
        count = 0
        with open(vcf_path, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    count += 1
        return count

    before = count_variants(background_vcf)

    # Compress and index the BED file for bcftools -T
    bed_gz = temp_dir / (Path(exclusion_bed).name + '.gz')
    subprocess.run(['bgzip', '-c', exclusion_bed],
                   stdout=open(bed_gz, 'wb'), check=True)
    subprocess.run(['tabix', '-p', 'bed', str(bed_gz)], check=True)

    # Compress and index the background VCF
    bg_gz = temp_dir / (Path(background_vcf).name + '.filter.gz')
    subprocess.run(['bgzip', '-c', background_vcf],
                   stdout=open(bg_gz, 'wb'), check=True)
    subprocess.run(['tabix', '-p', 'vcf', str(bg_gz)], check=True)

    # Exclude positions: bcftools view -T ^bed.gz
    cmd = [
        'bcftools', 'view',
        '-T', f'^{bed_gz}',   # ^ means exclude these regions
        '-O', 'v',             # Output uncompressed VCF
        '-o', output_vcf,
        str(bg_gz)
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"bcftools view -T (exclude) failed: {result.stderr}")

    after = count_variants(output_vcf)
    removed = before - after

    # Clean up
    bed_gz.unlink()
    (bed_gz.parent / (bed_gz.name + '.tbi')).unlink()
    bg_gz.unlink()
    (bg_gz.parent / (bg_gz.name + '.tbi')).unlink()

    return removed


def extract_sample_bcftools(background_vcf: str, sample_id: str, output_vcf: str):
    """Extract a single sample using bcftools"""
    cmd = [
        'bcftools', 'view',
        '-s', sample_id,
        '-O', 'v',  # Output uncompressed VCF
        '-o', output_vcf,
        background_vcf
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"bcftools view failed: {result.stderr}")


def merge_vcfs_bcftools(background_vcf: str, variant_vcf: str, output_vcf: str, temp_dir: Path):
    """
    Merge VCF files using bcftools concat.
    After pre-merge filtering, there should be no overlapping positions,
    so -a (allow overlaps) is kept as a safety net but should not trigger.
    """
    # First, compress and index both VCFs
    bg_gz = temp_dir / (Path(background_vcf).name + '.gz')
    var_gz = temp_dir / (Path(variant_vcf).name + '.gz')

    # Compress background VCF
    subprocess.run(['bgzip', '-c', background_vcf],
                   stdout=open(bg_gz, 'wb'), check=True)
    subprocess.run(['tabix', '-p', 'vcf', str(bg_gz)], check=True)

    # Compress variant VCF
    subprocess.run(['bgzip', '-c', variant_vcf],
                   stdout=open(var_gz, 'wb'), check=True)
    subprocess.run(['tabix', '-p', 'vcf', str(var_gz)], check=True)

    # Merge using bcftools concat (allows overlapping variants as safety net)
    cmd = [
        'bcftools', 'concat',
        '-a',  # Allow overlaps (safety net)
        '-O', 'v',  # Output uncompressed VCF
        '-o', output_vcf,
        str(bg_gz), str(var_gz)
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"bcftools concat failed: {result.stderr}")

    # Sort the merged VCF
    sorted_vcf = output_vcf + '.sorted'
    cmd = [
        'bcftools', 'sort',
        '-O', 'v',
        '-o', sorted_vcf,
        output_vcf
    ]
    subprocess.run(cmd, check=True)

    # Replace original with sorted
    os.replace(sorted_vcf, output_vcf)

    # Clean up temporary files
    bg_gz.unlink()
    (bg_gz.parent / (bg_gz.name + '.tbi')).unlink()
    var_gz.unlink()
    (var_gz.parent / (var_gz.name + '.tbi')).unlink()


def main():
    parser = argparse.ArgumentParser(
        description='Batch generate virtual patient VCF files using bcftools, '
                    'with deduplication and pathogenic-priority overlap handling. '
                    'When background and pathogenic VCFs share the same (CHROM, POS), '
                    'the background record is removed and the pathogenic record is kept.'
    )
    parser.add_argument(
        '--mapping',
        default='used_samples_combine.txt',
        help='Path to case-sample mapping file'
    )
    parser.add_argument(
        '--samples-dir',
        default='samples',
        help='Directory containing background sample VCF files'
    )
    parser.add_argument(
        '--variants-base-dir',
        default='variants_vcf_38',
        help='Base directory containing GRCh37 and GRCh38 subdirectories'
    )
    parser.add_argument(
        '--output-dir',
        default='genotype',
        help='Output directory for patient VCF files'
    )
    parser.add_argument(
        '--cases',
        nargs='+',
        help='Specific case numbers to process'
    )
    parser.add_argument(
        '--temp-dir',
        default='temp_vcfs',
        help='Temporary directory'
    )
    parser.add_argument(
        '--compress',
        action='store_true',
        help='Compress output VCF files with bgzip and create index'
    )
    parser.add_argument(
        '--skip-overlap-filter',
        action='store_true',
        help='Skip the pre-merge overlap filtering step. '
             'Without this flag, background sites at pathogenic positions '
             'are removed before merging (pathogenic takes priority).'
    )
    parser.add_argument(
        '--variants-data',
        default=None,
        help='Path to zygosity-annotated variants data file (.with_zygosity.txt). '
             'When provided, GT is set based on inheritance mode and zygosity '
             '(e.g., 1/1 for AR homozygous, 1 for X-linked hemizygous). '
             'Without this flag, all pathogenic variants default to GT=0/1.'
    )

    args = parser.parse_args()

    # Check for bcftools
    if not check_bcftools():
        print("ERROR: bcftools is not installed or not in PATH")
        print("Please install bcftools: https://samtools.github.io/bcftools/")
        sys.exit(1)

    # Create directories
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    temp_dir = Path(args.temp_dir)
    temp_dir.mkdir(exist_ok=True)

    # Read mapping
    print(f"Reading case-sample mapping from {args.mapping}...")
    mapping = read_sample_mapping(args.mapping)
    print(f"Found {len(mapping)} case-sample mappings")

    # Read zygosity data if provided
    zygosity_data = {}
    if args.variants_data:
        print(f"Reading zygosity data from {args.variants_data}...")
        zygosity_data = read_zygosity_data(args.variants_data)
        n_non_het = sum(1 for v in zygosity_data.values() if v['gt'] != '0/1')
        print(f"Found {len(zygosity_data)} cases with zygosity info "
              f"({n_non_het} with non-0/1 GT)")
    else:
        print("No --variants-data provided; all pathogenic variants will use GT=0/1")

    # Filter cases if specified
    if args.cases:
        mapping = {k: v for k, v in mapping.items() if k in args.cases}
        print(f"Processing {len(mapping)} specified cases")

    # Process each case
    success_count = 0
    error_count = 0
    total_variant_dedup = 0       # variant VCF 内部去重数
    total_bg_removed = 0          # 背景 VCF 中被排除的重叠位点数
    total_gt_override = 0         # GT 从默认 0/1 被修改的 case 数

    for case_num, sample_id in mapping.items():
        try:
            print(f"\nProcessing Case {case_num} (Sample: {sample_id})...")

            # Try to find variant VCF in GRCh37 or GRCh38 directory
            variant_vcf = None
            ref_genome = None

            for genome_version in ['GRCh37', 'GRCh38']:
                test_path = Path(args.variants_base_dir) / genome_version / f"case{case_num}.vcf"
                if test_path.exists():
                    variant_vcf = test_path
                    ref_genome = genome_version
                    print(f"  Found variant VCF in {genome_version} directory")
                    break

            if variant_vcf is None:
                print(f"  WARNING: Variant VCF not found in GRCh37 or GRCh38 directories")
                error_count += 1
                continue

            # Paths
            background_vcf = Path(args.samples_dir) / f"{sample_id}.merged.vcf.gz"

            temp_sample_vcf = temp_dir / f"temp_{sample_id}_case{case_num}.vcf"
            temp_variant_norm = temp_dir / f"temp_variant_case{case_num}.vcf"
            # Create genome version specific output directory
            genome_output_dir = output_dir / ref_genome
            genome_output_dir.mkdir(exist_ok=True)

            output_vcf = genome_output_dir / f"patient_case{case_num}_{sample_id}.vcf"

            # Check if output already exists
            if output_vcf.exists() or (genome_output_dir / (output_vcf.name + '.gz')).exists():
                print(f"  Skipping Case {case_num} (Sample: {sample_id}): Output file already exists in {ref_genome}")
                continue

            # Check if files exist
            if not background_vcf.exists():
                print(f"  WARNING: Background VCF not found: {background_vcf}")
                error_count += 1
                continue

            # Step 1: Extract sample using bcftools
            print(f"  Extracting sample {sample_id} using bcftools...")
            extract_sample_bcftools(str(background_vcf), sample_id, str(temp_sample_vcf))

            # Step 2: Normalize pathogenic variant VCF (with internal dedup)
            # Determine GT from zygosity data
            gt_value = "0/1"  # default
            if case_num in zygosity_data:
                gt_value = zygosity_data[case_num]['gt']
                inh_mode = zygosity_data[case_num]['inheritance_mode']
                zyg = zygosity_data[case_num]['zygosity']
                if gt_value != "0/1":
                    print(f"  [ZYGOSITY] {inh_mode} / {zyg} → GT={gt_value}")
                    total_gt_override += 1
            print(f"  Normalizing pathogenic variant VCF (with coordinate dedup, GT={gt_value})...")
            total_before, total_after = normalize_pathogenic_vcf(
                str(variant_vcf), str(temp_variant_norm), sample_id,
                gt_value=gt_value
            )
            variant_dups = total_before - total_after
            total_variant_dedup += variant_dups

            # Step 3: Pre-merge overlap filtering
            # Remove from background VCF any sites at pathogenic positions
            # so that pathogenic variant always takes priority
            if not args.skip_overlap_filter:
                print(f"  Filtering background VCF to exclude pathogenic positions...")

                # 3a: Extract pathogenic positions
                patho_positions = extract_pathogenic_positions(str(temp_variant_norm))
                print(f"    Pathogenic variant positions: {len(patho_positions)}")

                if patho_positions:
                    # 3b: Write exclusion BED file
                    exclusion_bed = temp_dir / f"exclude_case{case_num}.bed"
                    write_exclusion_bed(patho_positions, str(exclusion_bed))

                    # 3c: Filter background VCF — remove overlapping sites
                    temp_bg_filtered = temp_dir / f"temp_{sample_id}_case{case_num}_filtered.vcf"
                    bg_removed = filter_background_exclude_positions(
                        str(temp_sample_vcf), str(exclusion_bed),
                        str(temp_bg_filtered), temp_dir
                    )
                    if bg_removed > 0:
                        print(f"    [OVERLAP] Removed {bg_removed} background site(s) at "
                              f"pathogenic position(s) — pathogenic variant takes priority")
                        total_bg_removed += bg_removed
                    else:
                        print(f"    No overlapping sites found between background and pathogenic VCFs")

                    # Replace unfiltered background with filtered version
                    temp_sample_vcf.unlink()
                    temp_sample_vcf = temp_bg_filtered

                    # Clean up BED file
                    exclusion_bed.unlink()
                else:
                    print(f"    No pathogenic positions to exclude (empty variant VCF)")

            # Step 4: Merge using bcftools
            print(f"  Merging VCFs using bcftools...")
            merge_vcfs_bcftools(str(temp_sample_vcf), str(temp_variant_norm),
                              str(output_vcf), temp_dir)

            # Compress if requested
            if args.compress:
                print(f"  Compressing output VCF...")
                subprocess.run(['bgzip', '-f', str(output_vcf)], check=True)
                subprocess.run(['tabix', '-p', 'vcf', str(output_vcf) + '.gz'], check=True)
                print(f"  ✓ Successfully created: {output_vcf}.gz")
            else:
                print(f"  ✓ Successfully created: {output_vcf}")

            # Clean up temp files
            if temp_sample_vcf.exists():
                temp_sample_vcf.unlink()
            temp_variant_norm.unlink()

            success_count += 1

        except Exception as e:
            print(f"  ERROR processing Case {case_num}: {str(e)}")
            error_count += 1
            continue

    # Summary
    print(f"\n{'='*60}")
    print(f"Processing complete!")
    print(f"  Successful:              {success_count}")
    print(f"  Errors:                  {error_count}")
    print(f"  Variant dedup removed:   {total_variant_dedup}")
    print(f"  Background sites removed: {total_bg_removed}  (pathogenic takes priority)")
    print(f"  GT overrides (non-0/1):  {total_gt_override}  (inheritance-aware zygosity)")
    print(f"  Output directory:        {output_dir}")
    print(f"{'='*60}")

    # Clean up temp directory if empty
    try:
        temp_dir.rmdir()
    except:
        pass


if __name__ == '__main__':
    main()
