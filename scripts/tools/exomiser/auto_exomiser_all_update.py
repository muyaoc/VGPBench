import os
import glob
import subprocess

# --- CORE EXOMISER CONFIGURATION ---
EXOMISER_JAR = "exomiser-cli-14.0.0.jar"
MEMORY = "8g"         # Memory allocation
DATA_VERSION = "2406" # Your exact data version

# --- DIRECTORY PATHS ---
# Set GENETOOLS_BASE to your project root before running, e.g.:
#   export GENETOOLS_BASE=/path/to/genetools
_base = os.environ.get("GENETOOLS_BASE", "")
if not _base:
    raise EnvironmentError("Set GENETOOLS_BASE to your project root "
                           "(e.g. export GENETOOLS_BASE=/path/to/genetools)")
HPO_DIR = os.path.join(_base, "phenotype")
VCF_BASE_DIR = os.path.join(_base, "genotype")
RESULTS_BASE_DIR = os.path.join(_base, "exomiser_data", "results")
# -----------------------------------

def process_all_cases():
    # Calculate absolute path for the Exomiser data directory
    current_dir = os.path.abspath(os.path.dirname(__file__)).replace("\\", "/")
    absolute_data_dir = f"{current_dir}/data"

    if not os.path.exists(HPO_DIR):
        print(f"Error: HPO directory not found at {HPO_DIR}")
        return

    # Find all hpo files ending with _hpos.txt
    hpo_files = [f for f in os.listdir(HPO_DIR) if f.endswith("_hpos.txt")]
    
    if not hpo_files:
        print(f"No HPO files found in {HPO_DIR}")
        return
        
    print(f"[*] Found {len(hpo_files)} cases to process. Starting batch pipeline...\n")

    for hpo_file in hpo_files:
        # Extract the <case#>_<sample id> part (e.g., 'case1_HG00407')
        prefix = hpo_file.replace("_hpos.txt", "")
        hpo_path = os.path.join(HPO_DIR, hpo_file)
        
        print(f"==================================================")
        print(f"[*] Checking Case: {prefix}")

        # --- SETUP OUTPUT DIRECTORY & CHECK FOR EXISTING RESULTS ---
        output_dir = os.path.join(RESULTS_BASE_DIR, prefix).replace("\\", "/")
        output_file_name = f"{prefix}_exomiser"
        
        # Exomiser outputs TSV_GENE as <output_file_name>.genes.tsv
        expected_tsv_file = os.path.join(output_dir, f"{output_file_name}.genes.tsv")
        
        if os.path.exists(expected_tsv_file):
            print(f"[!] Skipping {prefix}: TSV_GENE results already exist.")
            continue
        
        os.makedirs(output_dir, exist_ok=True)
        # -----------------------------------------------------------
        
        # 1. Locate VCF & Detect Genome Assembly automatically
        vcf_grch37 = os.path.join(VCF_BASE_DIR, "GRCh37", f"patient_{prefix}.vcf.gz")
        vcf_grch38 = os.path.join(VCF_BASE_DIR, "GRCh38", f"patient_{prefix}.vcf.gz")
        
        if os.path.exists(vcf_grch37):
            vcf_path = vcf_grch37
            assembly = "hg19"
        elif os.path.exists(vcf_grch38):
            vcf_path = vcf_grch38
            assembly = "hg38"
        else:
            print(f"[-] Skipping {prefix}: Could not find matching VCF in GRCh37 or GRCh38 folders.")
            continue

        print(f"[*] Assembly Detected: {assembly}")
        
        # 2. Extract HPOs
        with open(hpo_path, 'r') as f:
            hpo_list =[line.strip() for line in f if line.strip().startswith("HP:")]
        
        if not hpo_list:
            print(f"[-] Skipping {prefix}: No HPO terms found.")
            continue
            
        hpo_formatted_string = "[" + ", ".join(f"'{hpo}'" for hpo in hpo_list) + "]"
        
        # 3. Format Paths for Java / YAML
        vcf_forward_slash = vcf_path.replace("\\", "/")
        
        # Save the YAML directly in the patient's output folder
        yaml_filename = os.path.join(output_dir, f"{prefix}_analysis.yml").replace("\\", "/")

        # 4. Generate the YAML mapping
        yaml_content = f"""---
analysis:
    genomeAssembly: {assembly}
    vcf: {vcf_forward_slash}
    hpoIds: {hpo_formatted_string}
    inheritanceModes: {{
      AUTOSOMAL_DOMINANT: 0.1,
      AUTOSOMAL_RECESSIVE_HOM_ALT: 0.1,
      AUTOSOMAL_RECESSIVE_COMP_HET: 2.0,
      X_DOMINANT: 0.1,
      X_RECESSIVE_HOM_ALT: 0.1,
      X_RECESSIVE_COMP_HET: 2.0,
      MITOCHONDRIAL: 0.2
    }}
    analysisMode: PASS_ONLY
    frequencySources: [ 
        UK10K, GNOMAD_E_AFR, GNOMAD_E_AMR, GNOMAD_E_EAS, GNOMAD_E_NFE, GNOMAD_E_SAS,
        GNOMAD_G_AFR, GNOMAD_G_AMR, GNOMAD_G_EAS, GNOMAD_G_NFE, GNOMAD_G_SAS
    ]
    pathogenicitySources: [ REVEL, MVP ]
    steps:[
        failedVariantFilter: {{}},
        variantEffectFilter: {{
          remove: [
              FIVE_PRIME_UTR_EXON_VARIANT, FIVE_PRIME_UTR_INTRON_VARIANT,
              THREE_PRIME_UTR_EXON_VARIANT, THREE_PRIME_UTR_INTRON_VARIANT,
              NON_CODING_TRANSCRIPT_EXON_VARIANT, NON_CODING_TRANSCRIPT_INTRON_VARIANT,
              CODING_TRANSCRIPT_INTRON_VARIANT, UPSTREAM_GENE_VARIANT,
              DOWNSTREAM_GENE_VARIANT, INTERGENIC_VARIANT, REGULATORY_REGION_VARIANT
          ]
        }},
        frequencyFilter: {{ maxFrequency: 2.0 }},
        pathogenicityFilter: {{ keepNonPathogenic: true }},
        inheritanceFilter: {{}},
        omimPrioritiser: {{}},
        hiPhivePrioritiser: {{}}
    ]

outputOptions:
    outputContributingVariantsOnly: false
    numGenes: 0
    outputDirectory: {output_dir}
    outputFileName: {output_file_name}
    outputFormats: [TSV_GENE]
"""
    # outputFormats: [HTML, JSON, TSV_GENE, TSV_VARIANT, VCF]
        with open(yaml_filename, 'w') as f:
            f.write(yaml_content)
            
        # 5. Launch Exomiser
        cmd =[
            "java", f"-Xmx{MEMORY}", "-jar", EXOMISER_JAR,
            "--analysis", yaml_filename,
            f"--exomiser.data-directory={absolute_data_dir}",
            f"--exomiser.hg19.data-version={DATA_VERSION}",
            f"--exomiser.hg38.data-version={DATA_VERSION}",
            f"--exomiser.phenotype.data-version={DATA_VERSION}"
        ]
        
        print(f"[*] Running Exomiser engine...")
        try:
            subprocess.run(cmd, check=True)
            print(f"[+] Success! Results saved in:\n    {output_dir}\n")
        except subprocess.CalledProcessError as e:
            print(f"[-] Exomiser encountered an error while processing {prefix}.\n")

if __name__ == "__main__":
    process_all_cases()