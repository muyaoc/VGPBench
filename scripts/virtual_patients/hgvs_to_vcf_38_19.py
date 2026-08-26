#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import json
import os
import re
import sys
import time
import subprocess
import tempfile
from collections import defaultdict
from typing import Dict, List, Tuple, Optional, Set

import requests

SERVER_38 = "https://rest.ensembl.org"
SERVER_37 = "https://grch37.rest.ensembl.org"

EXT_HGVS = "/vep/human/hgvs?vcf_string=1&refseq=1"

MAX_PER_POST_37 = 50
MAX_PER_POST_38 = 200

NM_VERSION_RE = re.compile(r"^(N[MR]_\d+)\.\d+:(.+)$")
C_WEIRD_SUB_RE = re.compile(r"^c\.([ACGT])(\d+)([ACGT])$")

HGVS_PREFIX_OK = re.compile(r"^[cgmnr]\.", re.IGNORECASE)

NON_HGVS_PATTERNS = [
    r"Exon", r"exon", r"and",
    r"_\-",
    r"\-_",
    r"^c\.\d+[_-]\d+$",
    r"->",
]

GENE_TYPO_FIX = {
    "GNOA1": "GNAO1",
    "TWSIT1": "TWIST1",
    "MLL": "KMT2A",
    "MLL2": "KMT2D",
    "KIAA2022": "NEXMIF",
}

GENE_PREFIX_RE = re.compile(r"^([A-Za-z0-9]+):(.+)$")

ERR_RUNNER_IN = "Bio::EnsEMBL::VEP::Runner::IN"
SPLIT_SEP_RE = re.compile(r"[,\uFF0C;\uFF1B]")
REF_MISMATCH_RE = re.compile(r"Reference allele extracted from .*? \(([ACGT]+)\) does not match reference allele given by HGVS notation .*? \(([ACGT]+)\)")


def parse_args():
    p = argparse.ArgumentParser(
        description="Convert HGVS to per-case VCF via Ensembl REST (resume + ENST fallbacks + robust splitting) with hg38->hg19 conversion."
    )
    p.add_argument("--merged", default="hg38_variants.txt")
    p.add_argument("--outdir", default="variants_vcf_38to19")
    p.add_argument("--outdir_hg19", default="variants_vcf_hg19", help="Output directory for hg19 VCF files")
    p.add_argument("--failures", default="hgvs_to_vcf_failures_new_20260325.txt")
    p.add_argument("--unconvertible", default="unconvertible_hgvs_new_20260325.txt")
    p.add_argument("--sleep", type=float, default=0.2)
    p.add_argument("--timeout", type=float, default=20.0)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--skip_hg19_conversion", action="store_true", help="Skip hg38 to hg19 conversion")
    p.add_argument("--liftover_path", default="./liftOver", help="Path to UCSC liftOver executable")
    p.add_argument("--chain_file", default="hg38ToHg19.over.chain.gz", help="Path to hg38ToHg19 chain file")
    
    p.add_argument("--enable_nm_to_enst", action="store_true", default=True)
    p.add_argument("--disable_nm_to_enst", action="store_true", default=False)
    p.add_argument("--enable_gene_to_enst", action="store_true", default=True)
    p.add_argument("--disable_gene_to_enst", action="store_true", default=False)
    p.add_argument("--max_enst_per_gene", type=int, default=1)
    p.add_argument("--case", help="Process only this specific case ID (e.g. '435')")
    p.add_argument("--case_list", help="File containing list of case IDs to process")
    return p.parse_args()


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def check_liftover_tools(args):
    """检查liftOver工具和chain文件是否存在"""
    if args.skip_hg19_conversion:
        return True, None
    
    # 检查liftOver可执行文件
    if not os.path.exists(args.liftover_path):
        # 尝试在PATH中查找
        import shutil
        liftover_path = shutil.which("liftOver")
        if liftover_path:
            args.liftover_path = liftover_path
        else:
            return False, f"liftOver executable not found at {args.liftover_path}. Please install UCSC liftOver or specify with --liftover_path"
    
    # 检查chain文件
    if not os.path.exists(args.chain_file):
        # 尝试从常见位置查找
        common_locations = [
            args.chain_file,
            "/usr/local/share/ucsc/hg38ToHg19.over.chain.gz",
            "/usr/share/ucsc/hg38ToHg19.over.chain.gz",
            "hg38ToHg19.over.chain.gz"
        ]
        for loc in common_locations:
            if os.path.exists(loc):
                args.chain_file = loc
                break
        else:
            # 如果不存在，尝试下载
            print(f"[WARNING] Chain file not found at {args.chain_file}")
            print("[INFO] You can download it with:")
            print("  wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz")
            return False, f"Chain file not found: {args.chain_file}"
    
    # 测试liftOver是否可以运行
    try:
        result = subprocess.run([args.liftover_path], capture_output=True, text=True)
        if "usage: liftOver" in result.stdout or "usage: liftOver" in result.stderr:
            return True, None
    except Exception as e:
        return False, f"liftOver test failed: {str(e)}"
    
    return True, None


def run_liftover(bed_records: List[Dict], liftover_path: str, chain_file: str) -> Tuple[List[Dict], List[Dict]]:
    """
    使用liftOver转换坐标
    输入: [{"chrom": "chr1", "pos": "123456", "ref": "A", "alt": "G", "hgvs": "NM_xxx:c.123A>G", ...}, ...]
    输出: (转换成功的记录列表, 转换失败的记录列表)
    """
    if not bed_records:
        return [], []
    
    # 创建临时文件
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as input_bed, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as output_bed, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as unmapped_bed:
        
        input_bed_path = input_bed.name
        output_bed_path = output_bed.name
        unmapped_bed_path = unmapped_bed.name
        
        # 写入BED文件 (0-based, half-open)
        for i, rec in enumerate(bed_records):
            try:
                chrom = rec["chrom"]
                if not chrom.startswith("chr"):
                    chrom = f"chr{chrom}"
                pos = int(rec["pos"]) - 1  # BED是0-based, VCF是1-based
                end = pos + len(rec["ref"])
                # 添加索引以便后续匹配
                input_bed.write(f"{chrom}\t{pos}\t{end}\t{i}\n")
            except (ValueError, KeyError) as e:
                print(f"[WARNING] Invalid record for liftOver: {rec}, error: {e}")
                continue
        
        input_bed.flush()
        
        # 运行liftOver
        cmd = [
            liftover_path,
            input_bed_path,
            chain_file,
            output_bed_path,
            unmapped_bed_path
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # 读取转换结果
            converted = {}
            if os.path.exists(output_bed_path) and os.path.getsize(output_bed_path) > 0:
                with open(output_bed_path, 'r') as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 4:
                            chrom_hg19 = parts[0]
                            start_hg19 = int(parts[1])  # 0-based
                            original_idx = int(parts[3])
                            # 转换回1-based坐标
                            pos_hg19 = start_hg19 + 1
                            converted[original_idx] = (chrom_hg19, pos_hg19)
            
            # 处理转换结果
            successful = []
            failed = []
            
            for i, rec in enumerate(bed_records):
                if i in converted:
                    chrom_hg19, pos_hg19 = converted[i]
                    # 移除可能的"chr"前缀以保持一致性
                    if chrom_hg19.startswith("chr"):
                        chrom_hg19 = chrom_hg19[3:]
                    
                    # 创建新的记录
                    new_rec = rec.copy()
                    new_rec["chrom"] = chrom_hg19
                    new_rec["pos"] = str(pos_hg19)
                    new_rec["original_hg38_chrom"] = rec["chrom"]
                    new_rec["original_hg38_pos"] = rec["pos"]
                    successful.append(new_rec)
                else:
                    failed.append(rec)
            
            return successful, failed
            
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] liftOver failed: {e}")
            print(f"[ERROR] stderr: {e.stderr}")
            return [], bed_records
        finally:
            # 清理临时文件
            for path in [input_bed_path, output_bed_path, unmapped_bed_path]:
                try:
                    if os.path.exists(path):
                        os.unlink(path)
                except:
                    pass


def convert_hg38_to_hg19(records_hg38: List[Dict], args) -> Tuple[List[Dict], List[Dict]]:
    """
    将hg38记录转换为hg19记录
    """
    if not records_hg38:
        return [], []
    
    print(f"[INFO] Converting {len(records_hg38)} variants from hg38 to hg19...")
    
    # 运行liftOver
    successful, failed = run_liftover(records_hg38, args.liftover_path, args.chain_file)
    
    conversion_rate = len(successful) / len(records_hg38) * 100 if records_hg38 else 0
    print(f"[INFO] Conversion successful: {len(successful)}/{len(records_hg38)} ({conversion_rate:.1f}%)")
    
    if failed:
        print(f"[WARNING] {len(failed)} variants failed to convert:")
        for rec in failed[:5]:  # 只显示前5个失败的
            print(f"  - {rec.get('hgvs', 'Unknown')}: {rec.get('chrom')}:{rec.get('pos')}")
        if len(failed) > 5:
            print(f"  ... and {len(failed) - 5} more")
    
    return successful, failed


def write_case_vcf(out_path: str, genome: str, records: List[Dict], is_hg19: bool = False):
    """
    写入VCF文件
    is_hg19: 是否为hg19版本，决定INFO字段的内容
    """
    ensure_dir(os.path.dirname(out_path))
    
    with open(out_path, "w", encoding="utf-8") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##source=EnsemblREST_VEP_HGVS\n")
        out.write(f"##reference_genome={genome}\n")
        
        if is_hg19:
            out.write("##INFO=<ID=HGVS,Number=1,Type=String,Description=\"Input HGVS notation used for Ensembl REST\">\n")
            out.write("##INFO=<ID=Original_hg38,Number=1,Type=String,Description=\"Original hg38 coordinates (chrom:pos)\">\n")
            info_field = "HGVS={hgvs};Original_hg38={original_hg38}"
        else:
            out.write("##INFO=<ID=HGVS,Number=1,Type=String,Description=\"Input HGVS notation used for Ensembl REST\">\n")
            info_field = "HGVS={hgvs}"
        
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        
        for rec in records:
            hgvs = (rec.get("hgvs") or "").replace(" ", "_")
            
            if is_hg19:
                original_hg38 = f"{rec.get('original_hg38_chrom', '')}:{rec.get('original_hg38_pos', '')}"
                info = info_field.format(hgvs=hgvs, original_hg38=original_hg38)
            else:
                info = info_field.format(hgvs=hgvs)
            
            out.write(f'{rec["chrom"]}\t{rec["pos"]}\t.\t{rec["ref"]}\t{rec["alt"]}\t.\t.\t{info}\n')


def case_vcf_path(outdir: str, genome: str, case: str) -> str:
    """生成VCF文件路径"""
    if genome == "GRCh19":
        genome_dir = "GRCh19"
    elif genome == "hg19":
        genome_dir = "hg19"
    else:
        genome_dir = genome
    return os.path.join(outdir, genome_dir, f"case{case}.vcf")


# ... [保持原有的辅助函数不变，直到main函数] ...

# 确保这些函数都在新脚本中定义

# ... 上面是新添加的liftOver相关函数 ...

# 下面是需要从原脚本复制的函数

def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def pick_server_for_genome(genome: str) -> Optional[str]:
    if genome == "GRCh37":
        return SERVER_37
    if genome == "GRCh38":
        return SERVER_38
    return None


def max_per_post_for_genome(genome: str) -> int:
    return MAX_PER_POST_37 if genome == "GRCh37" else MAX_PER_POST_38


def chunked(items: List[str], size: int):
    for i in range(0, len(items), size):
        yield items[i:i + size]


def normalize_refseq_id(refseq: str) -> str:
    r = (refseq or "").strip()
    if not r:
        return r
    r = re.sub(r"^(N[MR])(\d+)(\.\d+)?$", r"\1_\2\3", r)
    return r


def normalize_change(chgvs: str) -> str:
    s = (chgvs or "").strip()
    s = s.replace("–", "-").replace("—", "-")
    s = s.replace("；", ";").replace("，", ",")
    s = re.sub(r"\s+", "", s)
    s = re.sub(r"de([ACGT]+)$", r"del\1", s)
    m = C_WEIRD_SUB_RE.match(s)
    if m:
        s = f"c.{m.group(2)}{m.group(1)}>{m.group(3)}"
    s = re.sub(r"^c\.([ACGT])(\d+)([ACGT]>[ACGT]+)$", r"c.\2\1>\3", s)
    s = re.sub(r"^c\.([ACGT])(\d+)([ACGT])$", r"c.\2\1>\3", s)
    if re.match(r"^\d+(del|ins|dup|[ACGT])", s):
        s = f"c.{s}"
    return s


def is_obviously_unparseable_change(change: str) -> bool:
    if not change or change == "-":
        return True
    if not HGVS_PREFIX_OK.match(change):
        return True
    for pat in NON_HGVS_PATTERNS:
        if re.search(pat, change):
            return True
    return False


def split_multi_field(s: str) -> List[str]:
    if not s:
        return []
    out = []
    for x in SPLIT_SEP_RE.split(s):
        x = x.strip()
        if x:
            out.append(x)
    return out


def expand_refseq_candidates(refseq_field: str) -> List[str]:
    parts = []
    for x in split_multi_field(refseq_field):
        x = normalize_refseq_id(x).strip()
        if not x or x == "-":
            continue
        parts.append(x)
    return parts


def append_unconvertible(path: str, case: str, genome: str, gene: str, refseq: str, change: str, reason: str):
    with open(path, "a", encoding="utf-8") as out:
        out.write(f"{case}\t{genome}\t{gene}\t{refseq}\t{change}\t{reason}\n")


def build_hgvs_candidates(refseq_field: str, gene: str, change_raw: str) -> List[str]:
    gene = (gene or "").strip()
    gene = GENE_TYPO_FIX.get(gene, gene)
    change = normalize_change(change_raw)

    if not change or change == "-":
        return []

    if change.startswith("m."):
        return [f"MT:{change}"]

    if is_obviously_unparseable_change(change):
        return []

    cands: List[str] = []

    refseqs = expand_refseq_candidates(refseq_field)
    for ref in refseqs:
        cands.append(f"{ref}:{change}")
        if "." in ref:
            base = ref.split(".", 1)[0]
            if base and base != ref:
                cands.append(f"{base}:{change}")

    if (not refseqs) and gene and gene != "-":
        cands.append(f"{gene}:{change}")

    cands = [c for c in cands if not c.startswith("-:")]

    seen = set()
    uniq = []
    for c in cands:
        if c in seen:
            continue
        seen.add(c)
        uniq.append(c)
    return uniq


def post_json_with_retry(url: str, payload: dict, timeout: float, max_tries: int = 6):
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    last = None
    for i in range(max_tries):
        r = requests.post(url, headers=headers, data=json.dumps(payload), timeout=timeout)
        last = r
        if r.ok:
            return r
        txt = r.text or ""
        if ERR_RUNNER_IN in txt:
            time.sleep(2.0 * (2 ** i))
            continue
        break
    if last is not None:
        raise requests.HTTPError(f"{last.status_code} {last.text}", response=last)
    raise requests.HTTPError("POST failed without response")


def get_json(url: str, timeout: float):
    headers = {"Accept": "application/json"}
    r = requests.get(url, headers=headers, timeout=timeout)
    if not r.ok:
        raise requests.HTTPError(f"{r.status_code} {r.text}", response=r)
    return r.json()


def query_vep_hgvs(server: str, hgvs_list: List[str], timeout: float):
    url = server + EXT_HGVS
    payload = {"hgvs_notations": hgvs_list}
    r = post_json_with_retry(url, payload, timeout=timeout, max_tries=6)
    return r.json()


def collect_vcf_records_and_errors(data: List[dict]) -> Tuple[List[dict], List[Tuple[str, str]]]:
    recs: List[dict] = []
    errs: List[Tuple[str, str]] = []
    for item in data:
        v_input = item.get("input") or ""
        if "error" in item and item.get("error"):
            errs.append((v_input, str(item.get("error"))[:400]))
            continue
        vcf_string = item.get("vcf_string")
        if not vcf_string:
            errs.append((v_input, "missing vcf_string"))
            continue
        parts = vcf_string.split("-")
        if len(parts) != 4:
            errs.append((v_input, f"bad vcf_string format: {vcf_string}"))
            continue
        chrom, pos, ref, alt = parts
        recs.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "hgvs": v_input})
    return recs, errs


def robust_query_hgvs_list(server: str, genome: str, hgvs_list: List[str], timeout: float, sleep_s: float):
    records: List[dict] = []
    failed: List[Tuple[str, str]] = []

    def _run(lst: List[str]):
        if not lst:
            return
        try:
            data = query_vep_hgvs(server, lst, timeout=timeout)
            
            if not isinstance(data, list):
                raise ValueError("API returned non-list data")
            
            if len(data) < len(lst) and len(lst) > 1:
                mid = len(lst) // 2
                _run(lst[:mid])
                _run(lst[mid:])
                return

            returned_inputs = {item.get("input") for item in data if item.get("input")}
            for h in lst:
                if h not in returned_inputs:
                    failed.append((h, "silent failure: missing from API response"))

            recs, errs = collect_vcf_records_and_errors(data)
            records.extend(recs)
            for h, reason in errs:
                failed.append((h, reason))
            time.sleep(sleep_s)
            return
        except requests.HTTPError as e:
            msg = str(e)
            if len(lst) > 1:
                mid = len(lst) // 2
                _run(lst[:mid])
                _run(lst[mid:])
                return
            failed.append((lst[0], f"HTTPError: {msg[:400]}"))
            return
        except requests.RequestException as e:
            if len(lst) > 1:
                mid = len(lst) // 2
                _run(lst[:mid])
                _run(lst[mid:])
                return
            failed.append((lst[0], f"RequestException: {str(e)[:400]}"))
            return

    for chunk in chunked(hgvs_list, max_per_post_for_genome(genome)):
        _run(chunk)

    return records, failed


def fetch_gene_transcripts(server: str, symbol: str, timeout: float, cache: dict) -> List[dict]:
    key = (server, symbol)
    if key in cache:
        return cache[key]
    url = server + EXT_LOOKUP_SYMBOL.format(symbol=symbol)
    try:
        data = get_json(url, timeout=timeout)
    except requests.HTTPError:
        cache[key] = []
        return []
    tx = data.get("Transcript") or []
    cache[key] = tx if isinstance(tx, list) else []
    return cache[key]


def gene_to_enst_list(server: str, gene_symbol: str, timeout: float, cache: dict, max_enst: int) -> List[str]:
    tx = fetch_gene_transcripts(server, gene_symbol, timeout=timeout, cache=cache)
    if not tx:
        return []
    canonical, others = [], []
    for t in tx:
        tid = t.get("id")
        if not tid or not str(tid).startswith("ENST"):
            continue
        if t.get("is_canonical") in (1, True, "1", "true", "True"):
            canonical.append(tid)
        else:
            others.append(tid)
    seen, out = set(), []
    for tid in canonical + others:
        if tid in seen:
            continue
        seen.add(tid)
        out.append(tid)
        if len(out) >= max_enst:
            break
    return out


def xrefs_id(server: str, xid: str, timeout: float, cache: dict) -> List[dict]:
    key = (server, xid)
    if key in cache:
        return cache[key]
    url = server + EXT_XREFS_ID.format(xid=xid)
    try:
        data = get_json(url, timeout=timeout)
    except requests.HTTPError:
        cache[key] = []
        return []
    cache[key] = data if isinstance(data, list) else []
    return cache[key]


def nm_to_enst_candidates(server: str, nm: str, timeout: float, cache: dict) -> List[str]:
    nm = normalize_refseq_id(nm)
    ids_to_try = [nm]
    if "." in nm:
        ids_to_try.append(nm.split(".", 1)[0])

    out, seen = [], set()
    for xid in ids_to_try:
        for rec in xrefs_id(server, xid, timeout=timeout, cache=cache):
            enst = rec.get("id")
            if not enst or not str(enst).startswith("ENST"):
                continue
            if enst in seen:
                continue
            seen.add(enst)
            out.append(enst)
    return out


def retry_with_nm_to_enst(server: str, genome: str, hgvs_list: List[str], timeout: float, sleep_s: float, xref_cache: dict):
    transformed, changed = [], 0
    for h in hgvs_list:
        if ":" not in h:
            transformed.append(h)
            continue
        prefix, change = h.split(":", 1)
        if prefix.startswith(("NM", "NR")):
            ensts = nm_to_enst_candidates(server, prefix, timeout=timeout, cache=xref_cache)
            if ensts:
                for enst in ensts:
                    transformed.append(f"{enst}:{change}")
                changed += 1
            else:
                transformed.append(h)
        else:
            transformed.append(h)

    seen, uniq = set(), []
    for t in transformed:
        if t in seen:
            continue
        seen.add(t)
        uniq.append(t)

    rec, bad = robust_query_hgvs_list(server, genome, uniq, timeout=timeout, sleep_s=sleep_s)
    
    if not rec and not bad:
        return [], [(h, "no NM/NR translation found") for h in hgvs_list]

    return rec, bad


def retry_with_gene_to_enst(server: str, genome: str, hgvs_list: List[str], timeout: float, sleep_s: float, gene_cache: dict, max_enst_per_gene: int):
    expanded, expanded_any = [], 0
    for h in hgvs_list:
        m = GENE_PREFIX_RE.match(h)
        if not m:
            expanded.append(h)
            continue
        prefix, change = m.group(1), m.group(2)
        if prefix.startswith(("NM", "NR", "ENST", "MT")):
            expanded.append(h)
            continue
        gene_symbol = GENE_TYPO_FIX.get(prefix, prefix)
        ensts = gene_to_enst_list(server, gene_symbol, timeout=timeout, cache=gene_cache, max_enst=max_enst_per_gene)
        if ensts:
            for enst in ensts:
                expanded.append(f"{enst}:{change}")
            expanded_any += 1
        else:
            expanded.append(h)

    seen, uniq = set(), []
    for t in expanded:
        if t in seen:
            continue
        seen.add(t)
        uniq.append(t)

    rec, bad = robust_query_hgvs_list(server, genome, uniq, timeout=timeout, sleep_s=sleep_s)
    if not rec and not bad:
        return [], [(h, "no ENST translation found for gene") for h in hgvs_list]
    
    return rec, bad


def retry_with_ref_mismatch_fix(server: str, genome: str, failed_items: List[Tuple[str, str]], timeout: float, sleep_s: float):
    candidates = []
    
    for hgvs, reason in failed_items:
        m = REF_MISMATCH_RE.search(reason)
        if not m:
            continue
        
        actual_ref = m.group(1)
        claimed_ref = m.group(2)
        
        sub_m = re.search(r"c\.(\d+)([ACGT])>([ACGT])", hgvs)
        if not sub_m:
            continue
            
        pos_str = sub_m.group(1)
        old_ref = sub_m.group(2)
        old_alt = sub_m.group(3)
        
        if actual_ref != old_alt:
            new_hgvs_1 = hgvs.replace(f"c.{pos_str}{old_ref}>{old_alt}", f"c.{pos_str}{actual_ref}>{old_alt}")
            candidates.append(new_hgvs_1)
            
        if actual_ref == old_alt:
            new_hgvs_2 = hgvs.replace(f"c.{pos_str}{old_ref}>{old_alt}", f"c.{pos_str}{actual_ref}>{old_ref}")
            candidates.append(new_hgvs_2)

    if not candidates:
        return [], []

    candidates = list(set(candidates))
    print(f"[INFO] Attempting Ref Mismatch Fix for {len(candidates)} candidates...")
    
    return robust_query_hgvs_list(server, genome, candidates, timeout=timeout, sleep_s=sleep_s)


def load_hgvs_by_case(merged_path: str, unconvertible_path: str) -> Dict[Tuple[str, str], List[str]]:
    if not os.path.isfile(merged_path):
        sys.stderr.write(f"[ERROR] merged file not found: {merged_path}\n")
        sys.exit(1)

    by_case: Dict[Tuple[str, str], List[str]] = defaultdict(list)

    with open(merged_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            case = str(row.get("#Case", "")).strip()
            genome = (row.get("Reference_genome") or "").strip()
            gene = (row.get("Gene") or "").strip()
            refseq_field = (row.get("Ref seq") or "").strip()
            chgvs_raw = (row.get("cHGVS") or "").strip()

            if not case or not genome:
                continue
            if not chgvs_raw or chgvs_raw == "-":
                continue

            gene_fixed = GENE_TYPO_FIX.get(gene, gene)

            for c in split_multi_field(chgvs_raw):
                if not c or c == "-":
                    continue

                change_norm = normalize_change(c)
                if change_norm and not change_norm.startswith("m.") and is_obviously_unparseable_change(change_norm):
                    append_unconvertible(
                        unconvertible_path, case, genome, gene_fixed, refseq_field, c,
                        "filtered_nonstandard_hgvs"
                    )
                    continue

                cands = build_hgvs_candidates(refseq_field, gene_fixed, c)
                if not cands:
                    append_unconvertible(
                        unconvertible_path, case, genome, gene_fixed, refseq_field, c,
                        "no_candidates_after_cleaning"
                    )
                    continue

                by_case[(case, genome)].extend(cands)

    for k in list(by_case.keys()):
        seen, uniq = set(), []
        for h in by_case[k]:
            if h in seen:
                continue
            seen.add(h)
            uniq.append(h)
        by_case[k] = uniq

    return by_case


def is_nonempty_file(path: str) -> bool:
    return os.path.isfile(path) and os.path.getsize(path) > 0


def append_failure(fail_path: str, case: str, genome: str, hgvs_count: int, reason: str):
    with open(fail_path, "a", encoding="utf-8") as out:
        out.write(f"{case}\t{genome}\t{hgvs_count}\t{reason}\n")


def main():
    args = parse_args()
    
    # 检查liftOver工具
    if not args.skip_hg19_conversion:
        liftover_ok, error_msg = check_liftover_tools(args)
        if not liftover_ok:
            print(f"[ERROR] {error_msg}")
            print("[INFO] You can skip hg19 conversion with --skip_hg19_conversion")
            sys.exit(1)
    
    enable_nm_to_enst = args.enable_nm_to_enst and (not args.disable_nm_to_enst)
    enable_gene_to_enst = args.enable_gene_to_enst and (not args.disable_gene_to_enst)

    merged_path = os.path.abspath(args.merged)
    outdir = os.path.abspath(args.outdir)
    outdir_hg19 = os.path.abspath(args.outdir_hg19) if not args.skip_hg19_conversion else None
    fail_path = os.path.abspath(args.failures)
    unconvertible_path = os.path.abspath(args.unconvertible)

    ensure_dir(outdir)
    if outdir_hg19:
        ensure_dir(outdir_hg19)
        
    if not os.path.exists(unconvertible_path) or os.path.getsize(unconvertible_path) == 0:
        with open(unconvertible_path, "w", encoding="utf-8") as out:
            out.write("Case\tGenome\tGene\tRefSeq\tcHGVS\tReason\n")

    print(f"[INFO] Merged table : {merged_path}")
    print(f"[INFO] Output dir (hg38): {outdir}")
    if not args.skip_hg19_conversion:
        print(f"[INFO] Output dir (hg19): {outdir_hg19}")
        print(f"[INFO] liftOver path: {args.liftover_path}")
        print(f"[INFO] Chain file: {args.chain_file}")
    else:
        print(f"[INFO] hg19 conversion: SKIPPED")
    print(f"[INFO] Fail log     : {fail_path}")
    print(f"[INFO] Unconvertible: {unconvertible_path}")
    print(f"[INFO] Overwrite    : {args.overwrite}")
    print(f"[INFO] NM->ENST      : {enable_nm_to_enst}")
    print(f"[INFO] Gene->ENST    : {enable_gene_to_enst} (max {args.max_enst_per_gene}/gene)")
    print()

    by_case = load_hgvs_by_case(merged_path, unconvertible_path)

    # Filter by --case if provided
    if args.case:
        target_case = str(args.case).strip()
        print(f"[INFO] Filtering for Case: {target_case}")
        filtered_cases = {}
        for k, v in by_case.items():
            if k[0] == target_case:
                filtered_cases[k] = v
        by_case = filtered_cases
        if not by_case:
            print(f"[WARN] Case {target_case} not found in loaded data.")
            return

    # Filter by --case_list if provided
    if args.case_list:
        list_path = os.path.abspath(args.case_list)
        if not os.path.exists(list_path):
            sys.stderr.write(f"[ERROR] Case list file not found: {list_path}\n")
            sys.exit(1)
        
        print(f"[INFO] Filtering using case list: {list_path}")
        allowed_cases = set()
        with open(list_path, 'r', encoding='utf-8') as f:
            for line in f:
                c = line.strip()
                if c:
                    allowed_cases.add(c)
        
        filtered_cases = {}
        found_count = 0
        for k, v in by_case.items():
            if k[0] in allowed_cases:
                filtered_cases[k] = v
                found_count += 1
        
        by_case = filtered_cases
        print(f"[INFO] Found {found_count} matching cases from list (out of {len(allowed_cases)} in list).")

    total_groups = len(by_case)

    skipped = 0
    success = 0
    failed = 0
    hg19_success = 0
    hg19_failed = 0

    gene_lookup_cache = {}
    xref_cache = {}

    def _case_key(x):
        (case, genome) = x[0]
        try:
            return (int(case), genome)
        except Exception:
            return (case, genome)

    for (case, genome), hgvs_list in sorted(by_case.items(), key=_case_key):
        out_path_hg38 = case_vcf_path(outdir, genome, case)
        out_path_hg19 = case_vcf_path(outdir_hg19, "GRCh19", case) if outdir_hg19 else None

        # 检查是否跳过
        skip_hg38 = (not args.overwrite) and is_nonempty_file(out_path_hg38)
        skip_hg19 = outdir_hg19 and (not args.overwrite) and is_nonempty_file(out_path_hg19)
        
        if skip_hg38 and (skip_hg19 or args.skip_hg19_conversion):
            skipped += 1
            continue

        server = pick_server_for_genome(genome)
        if server is None:
            failed += 1
            append_failure(fail_path, case, genome, len(hgvs_list), f"unsupported Reference_genome: {genome}")
            continue

        print(f"=== Case {case} / {genome} ===")
        print(f"[INFO] HGVS count: {len(hgvs_list)}")
        if not skip_hg38:
            print(f"[INFO] Output (hg38): {out_path_hg38}")
        if outdir_hg19 and not skip_hg19:
            print(f"[INFO] Output (hg19): {out_path_hg19}")

        # 1. 获取hg38坐标
        all_records, current_bad = robust_query_hgvs_list(
            server, genome, hgvs_list, timeout=args.timeout, sleep_s=args.sleep
        )

        # 2. 回退策略
        if current_bad and enable_nm_to_enst:
            to_retry = [h for h, r in current_bad]
            print(f"[INFO] Fallback (NM->ENST) for {len(to_retry)} items...")
            rec2, bad2 = retry_with_nm_to_enst(
                server, genome, to_retry, timeout=args.timeout, sleep_s=args.sleep, xref_cache=xref_cache
            )
            all_records.extend(rec2)
            current_bad = bad2

        if current_bad and enable_gene_to_enst:
            to_retry = [h for h, r in current_bad]
            print(f"[INFO] Fallback (Gene->ENST) for {len(to_retry)} items...")
            rec3, bad3 = retry_with_gene_to_enst(
                server, genome, to_retry, timeout=args.timeout, sleep_s=args.sleep,
                gene_cache=gene_lookup_cache, max_enst_per_gene=args.max_enst_per_gene
            )
            all_records.extend(rec3)
            current_bad = bad3

        if current_bad:
            has_mismatch = any(REF_MISMATCH_RE.search(r) for h, r in current_bad)
            if has_mismatch:
                print(f"[INFO] Fallback (RefMismatchFix) for {len(current_bad)} items...")
                rec4, bad4 = retry_with_ref_mismatch_fix(
                    server, genome, current_bad, timeout=args.timeout, sleep_s=args.sleep
                )
                all_records.extend(rec4)
                current_bad = bad4

        # 记录失败
        for h, reason in current_bad:
            append_failure(fail_path, case, genome, 1, f"{reason} | HGVS={h}")

        if not all_records:
            failed += 1
            print("[FAIL] no valid records (missing vcf_string)")
            print()
            continue

        # 去重
        uniq_key = set()
        uniq_records = []
        for r in all_records:
            k = (r["chrom"], r["pos"], r["ref"], r["alt"])
            if k in uniq_key:
                continue
            uniq_key.add(k)
            uniq_records.append(r)

        # 写入hg38 VCF
        if not skip_hg38:
            write_case_vcf(out_path_hg38, genome, uniq_records, is_hg19=False)
            print(f"[OK] wrote {len(uniq_records)} unique variants (hg38)")

        # 转换为hg19
        if not args.skip_hg19_conversion and genome.upper() in ["GRCH38", "HG38"]:
            hg19_records, hg19_failed_records = convert_hg38_to_hg19(uniq_records, args)
            
            if hg19_records and not skip_hg19:
                write_case_vcf(out_path_hg19, "GRCh19", hg19_records, is_hg19=True)
                hg19_success += 1
                print(f"[OK] wrote {len(hg19_records)} variants (hg19)")
            
            if hg19_failed_records:
                hg19_failed += 1
                print(f"[WARNING] {len(hg19_failed_records)} variants failed hg19 conversion")
        elif not args.skip_hg19_conversion:
            print(f"[INFO] Skipping hg19 conversion for genome: {genome}")

        success += 1
        print()

    print("========== SUMMARY ==========")
    print(f"Total (Case,Genome) groups : {total_groups}")
    print(f"Skipped (already exists)    : {skipped}")
    print(f"Success (hg38 generated)    : {success}")
    print(f"Failed                      : {failed}")
    if not args.skip_hg19_conversion:
        print(f"hg19 conversion successful : {hg19_success}")
        print(f"hg19 conversion failed     : {hg19_failed}")
    print(f"Failure log                 : {fail_path}")
    print(f"Unconvertible log           : {unconvertible_path}")
    print("=============================")


if __name__ == "__main__":
    main()