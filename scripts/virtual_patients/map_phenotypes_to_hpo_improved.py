#!/usr/bin/env python3
"""
Improved version: Map phenotype terms to HPO IDs with better error handling.
Uses local HPO mapping when possible and has better retry logic.
"""

import os
import re
import json
import time
from pathlib import Path
import urllib.request
import urllib.parse
import urllib.error
import socket
import argparse

# Common phenotype to HPO mappings (pre-populated to reduce API calls)
COMMON_HPO_MAPPINGS = {
    "Abnormality of the eye": "HP:0000478",
    "Abnormality of the nervous system": "HP:0000707",
    "Abnormality of the musculoskeletal system": "HP:0000924",
    "Multiple congenital anomalies": "HP:0001197",
    "Abnormality of the cardiovascular system": "HP:0001626",
    "Abnormality of the ear": "HP:0000598",
    "Abnormality of the skeletal system": "HP:0000924",
    "Aplasia/hypoplasia of the extremities": "HP:0009825",
    "Abnormality of the integument": "HP:0001574",
    "Abnormality of connective tissue": "HP:0003549",
    "Abnormality of the genitourinary system": "HP:0000119",
    "Abnormality of the endocrine system": "HP:0000818",
    "Abnormality of metabolism/homeostasis": "HP:0001939",
    "Abnormality of the immune system": "HP:0002715",
    "Abnormality of blood and blood-forming tissues": "HP:0001871",
    "Abnormality of the abdomen": "HP:0001438",
    "Abnormality of the mitochondrion": "HP:0001427",
    "Abnormality of the musculature": "HP:0003011",
    "Abnormality of the peripheral nervous system": "HP:0000759",
    "Growth abnormality": "HP:0001507",
    "Seizures": "HP:0001250",
    "Autism spectrum disorder": "HP:0000729",
    "Autism": "HP:0000717",
    "Global developmental delay": "HP:0001263",
    "Intellectual disability": "HP:0001249",
    "Delayed speech and language development": "HP:0000750",
    "Motor delay": "HP:0001270",
    "Muscular hypotonia": "HP:0001252",
    "Microcephaly": "HP:0000252",
    "Macrocephaly": "HP:0000256",
    "Nystagmus": "HP:0000639",
    "Ataxia": "HP:0001251",
    "Spasticity": "HP:0001257",
    "Dystonia": "HP:0001332",
    "Chorea": "HP:0002072",
    "Short stature": "HP:0004322",
    "Failure to thrive": "HP:0001508",
    "Vomiting": "HP:0002013",
    "Constipation": "HP:0002019",
    "Diarrhea": "HP:0002014",
    "Hepatomegaly": "HP:0002240",
    "Splenomegaly": "HP:0001744",
    "Hearing impairment": "HP:0000365",
    "Visual impairment": "HP:0000505",
    "Muscle weakness": "HP:0001324",
    "Hyperreflexia": "HP:0001347",
    "Hypoglycemia": "HP:0001943",
    "Metabolic acidosis": "HP:0001942",
    "Lactic acidosis": "HP:0003128",
    "Developmental regression": "HP:0002376",
    "Febrile seizures": "HP:0002373",
    "Generalized seizures": "HP:0002197",
    "Hypertonia": "HP:0001276",
    "Cryptorchidism": "HP:0000028",
    "Joint hypermobility": "HP:0001382",
    "Scoliosis": "HP:0002650",
    "Optic atrophy": "HP:0000648",
    "Cerebellar hypoplasia": "HP:0001321",
    "Brain atrophy": "HP:0012444",
}

def extract_hpo_ids_from_text(phenotype_text):
    """Extract HPO IDs that are already in HP:XXXXXXX format."""
    if not phenotype_text or phenotype_text.strip() == '-':
        return []
    
    hpo_pattern = r'HP:\s*(\d{7})'
    matches = re.findall(hpo_pattern, phenotype_text)
    return sorted(set([f"HP:{match}" for match in matches]))

def query_hpo_api(phenotype_term, max_retries=3):
    """Query HPO API to find HPO ID for a phenotype term with improved error handling."""
    base_url = "https://ontology.jax.org/api/hp/search"
    
    for attempt in range(max_retries):
        try:
            # Set longer timeout
            socket.setdefaulttimeout(30)
            
            # URL encode the search term
            params = urllib.parse.urlencode({'q': phenotype_term, 'max': 1})
            url = f"{base_url}?{params}"
            
            # Make request with headers
            req = urllib.request.Request(url)
            req.add_header('Accept', 'application/json')
            req.add_header('User-Agent', 'Mozilla/5.0')
            
            with urllib.request.urlopen(req, timeout=30) as response:
                data = json.loads(response.read().decode())
                
                # Extract HPO ID from response
                if data and 'terms' in data and len(data['terms']) > 0:
                    hpo_id = data['terms'][0].get('id', '')
                    if hpo_id:
                        return hpo_id
            
            return None
            
        except urllib.error.HTTPError as e:
            if e.code == 429:  # Rate limit
                wait_time = min(30, 2 ** attempt)
                print(f"  Rate limited, waiting {wait_time}s...")
                time.sleep(wait_time)
            else:
                print(f"  HTTP Error {e.code}")
                return None
        except (urllib.error.URLError, socket.timeout, OSError) as e:
            if attempt < max_retries - 1:
                wait_time = min(10, 2 ** attempt)
                print(f"  Network error (attempt {attempt + 1}/{max_retries}), waiting {wait_time}s...")
                time.sleep(wait_time)
            else:
                print(f"  Failed after {max_retries} attempts")
                return None
        except Exception as e:
            print(f"  Unexpected error: {type(e).__name__}")
            return None
    
    return None

def map_phenotypes_to_hpo(phenotype_text, cache):
    """Map phenotype text to HPO IDs."""
    # First check if already in HPO format
    existing_hpo = extract_hpo_ids_from_text(phenotype_text)
    if existing_hpo:
        return existing_hpo, []
    
    if not phenotype_text or phenotype_text.strip() == '-':
        return [], []
    
    # Split by comma or semicolon to get individual phenotype terms
    phenotype_terms = [term.strip() for term in re.split(r'[;,]', phenotype_text)]
    
    hpo_ids = []
    failed_terms = []
    for term in phenotype_terms:
        if not term:
            continue
        
        # Check common mappings first
        if term in COMMON_HPO_MAPPINGS:
            hpo_id = COMMON_HPO_MAPPINGS[term]
            hpo_ids.append(hpo_id)
            print(f"  Found (local): {term[:60]} -> {hpo_id}")
            continue
        
        # Check cache
        if term in cache:
            if cache[term]:
                hpo_ids.append(cache[term])
                print(f"  Found (cache): {term[:60]} -> {cache[term]}")
            else:
                print(f"  Found (cache): {term[:60]} -> None (Previous lookup failed)")
                if not cache[term]:
                     failed_terms.append((term, "Cached Failure"))
            continue
        
        # Query API with delay
        print(f"  Querying API: {term[:60]}...")
        time.sleep(1.5)  # Increased delay between requests
        
        hpo_id = query_hpo_api(term)
        
        if hpo_id:
            cache[term] = hpo_id
            hpo_ids.append(hpo_id)
            print(f"    Found: {hpo_id}")
        else:
            cache[term] = None
            failed_terms.append((term, "Not Found / API Error"))
            print(f"    Not found")
    
    return sorted(set(hpo_ids)), failed_terms

def load_sample_mapping(new_file):
    """Load case to sample mapping from reference files."""
    case_to_sample = {}
    
    with open(new_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                case_to_sample[parts[0]] = parts[1]
    
    return case_to_sample

# def load_sample_mapping(combine_file, new_file):
#     """Load case to sample mapping from reference files."""
#     case_to_sample = {}
    
#     with open(combine_file, 'r', encoding='utf-8') as f:
#         for line in f:
#             line = line.strip()
#             if not line or line.startswith('Case'):
#                 continue
#             parts = line.split()
#             if len(parts) >= 2:
#                 case_to_sample[parts[0]] = parts[1]
    
#     with open(new_file, 'r', encoding='utf-8') as f:
#         for line in f:
#             line = line.strip()
#             if not line:
#                 continue
#             parts = line.split()
#             if len(parts) >= 2:
#                 case_to_sample[parts[0]] = parts[1]
    
#     return case_to_sample

def process_filtered_data(filepath, case_to_sample, output_dir, cache_file, failure_file, target_case=None, ignore_cache=False):
    """Process filtered_cHGVS_data.txt with phenotype mapping."""
    
    # Load cache if exists and not ignored
    cache = {}
    if not ignore_cache and os.path.exists(cache_file):
        print(f"Loading cache from {cache_file}...")
        with open(cache_file, 'r', encoding='utf-8') as f:
            cache = json.load(f)
        print(f"Loaded {len(cache)} cached mappings")
    elif ignore_cache:
        print("Ignoring cache as requested.")
    
    created_files = []
    skipped_cases = []
    
    # Initialize failure file
    with open(failure_file, 'w', encoding='utf-8') as f:
        f.write("#Case\tPhenotype_Term\tDisease\tError\n")

    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#Case'):
                continue
            
            parts = line.split('\t')
            if len(parts) < 6:
                continue
            
            case_num = parts[0]
            
            # If a specific case is requested, skip others
            if target_case and case_num != target_case:
                continue

            phenotype = parts[5] if len(parts) > 5 else ''
            disease = parts[6] if len(parts) > 6 else 'Unknown'

            print(f"\nProcessing case {case_num}...")
            
            # Map phenotypes to HPO IDs
            hpo_ids, failed_terms = map_phenotypes_to_hpo(phenotype, cache)
            
            # Log failures
            if failed_terms:
                with open(failure_file, 'a', encoding='utf-8') as ff:
                    for term, error in failed_terms:
                        ff.write(f"{case_num}\t{term}\t{disease}\t{error}\n")

            if not hpo_ids:
                print(f"  [SKIP LOG] Case {case_num}: No HPO IDs found for phenotype: '{phenotype}'")
                skipped_cases.append(f"{case_num} (No HPO IDs)")
                continue
            
            # Get sample ID
            sample_id = case_to_sample.get(case_num)
            if not sample_id:
                print(f"  [SKIP LOG] Case {case_num}: No sample ID found in mapping files")
                skipped_cases.append(f"{case_num} (No Sample ID)")
                continue
            
            # Create output filename
            output_filename = f"case{case_num}_{sample_id}_hpos.txt"
            output_path = os.path.join(output_dir, output_filename)
            
            # Write HPO IDs to file
            with open(output_path, 'w', encoding='utf-8') as out_f:
                for hpo_id in hpo_ids:
                    out_f.write(f"{hpo_id}\n")
            
            created_files.append(output_filename)
            print(f"  Created: {output_filename} ({len(hpo_ids)} HPO IDs)")
            
            # Save cache every 5 files
            if len(created_files) % 5 == 0:
                with open(cache_file, 'w', encoding='utf-8') as f:
                    json.dump(cache, f, indent=2)
                print(f"  [Cache saved: {len(cache)} entries]")
    
    # Save final cache
    with open(cache_file, 'w', encoding='utf-8') as f:
        json.dump(cache, f, indent=2)
    print(f"\nSaved cache to {cache_file}")
    
    return created_files, skipped_cases

def main():
    parser = argparse.ArgumentParser(description='Map phenotypes to HPO IDs.')
    parser.add_argument('--case', type=str, help='Process only a specific case number')
    parser.add_argument('--ignore-cache', action='store_true', help='Ignore existing cache and re-query API')
    args = parser.parse_args()

    base_dir = Path(__file__).parent
    # filtered_file = base_dir / "filtered_cHGVS_data.txt"

    # filtered_file = base_dir / "published_diagnosed_pathogenic_variants.merged.txt"
    # combine_file = base_dir / "used_samples_combine.txt"
    # new_file = base_dir / "used_samples_new.txt"
    # output_dir = base_dir / "phenotype"
    # cache_file = base_dir / "phenotype_hpo_cache.json"
    # failure_file = base_dir / "phenotype_mapping_failures.txt"

    filtered_file = base_dir / "filtered_pheno_cases_data.txt"
    new_file = base_dir / "used_samples_new_novcf.txt"
    output_dir = base_dir / "phenotype_novcf"
    cache_file = base_dir / "phenotype_hpo_cache.json"
    failure_file = base_dir / "phenotype_mapping_failures_novcf.txt"
    
    # Create output directory
    output_dir.mkdir(exist_ok=True)
    
    # Load case to sample mapping
    print("Loading sample mappings...")
    # case_to_sample = load_sample_mapping(combine_file, new_file)
    case_to_sample = load_sample_mapping(new_file)
    print(f"Loaded {len(case_to_sample)} case-to-sample mappings\n")
    
    # Process filtered data with API mapping
    print("Processing filtered_cHGVS_data.txt with HPO mapping...")
    print(f"Using {len(COMMON_HPO_MAPPINGS)} pre-loaded common mappings")
    print("=" * 60)
    files, skipped = process_filtered_data(filtered_file, case_to_sample, output_dir, cache_file, failure_file, args.case, args.ignore_cache)
    
    # Summary
    print(f"\n{'='*60}")
    print(f"Summary:")
    print(f"  Total HPO files created: {len(files)}")
    print(f"  Skipped cases (no HPO IDs): {len(skipped)}")
    print(f"  Output directory: {output_dir}")
    print(f"  Cache file: {cache_file}")
    print(f"  Failures log: {failure_file}")
    print(f"{'='*60}")
    
    if skipped:
        print(f"\nSkipped cases: {', '.join(skipped[:20])}")
        if len(skipped) > 20:
            print(f"  ... and {len(skipped) - 20} more")

if __name__ == "__main__":
    main()
