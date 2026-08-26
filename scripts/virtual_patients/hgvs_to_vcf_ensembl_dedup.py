#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
批量把 published_diagnosed_pathogenic_variants.merged.txt 里的 HGVS
通过 Ensembl REST API 转成每个 Case 的 VCF：

GRCh37 → pathogenic_variants/GRCh37/case<Case>.vcf
GRCh38 → pathogenic_variants/GRCh38/case<Case>.vcf

相比 hgvs_to_vcf_ensembl.py，本脚本增加：
  1. 按 (CHROM, POS, REF, ALT) 坐标去重，避免不同转录本 HGVS 映射到
     同一基因组坐标时产生重复行
  2. 去重时保留第一条记录的 HGVS 信息，其余以 INFO 字段追加

依赖：
  pip install requests

注意：
  - 需要能访问 Ensembl REST（外网）
  - Ensembl 每次 POST 最多 200 条 HGVS，所以脚本会自动分批
"""

import argparse
import csv
import json
import os
import sys
import time
from collections import defaultdict, OrderedDict

import requests

# Ensembl REST endpoint
SERVER_38 = "https://rest.ensembl.org"
SERVER_37 = "https://grch37.rest.ensembl.org"

EXT_HGVS = "/vep/human/hgvs?vcf_string=1"  # 要求输出里带 vcf_string 字段

MAX_PER_POST = 200  # 每次最多 200 条 JSON（官方限制）


def parse_args():
    p = argparse.ArgumentParser(
        description="Use Ensembl REST VEP (HGVS endpoint) to convert HGVS to VCF per Case. "
                    "Deduplicates records by (CHROM, POS, REF, ALT) coordinates."
    )
    p.add_argument(
        "--merged",
        default="published_diagnosed_pathogenic_variants.merged.txt",
        help="Path to merged pathogenic table (TSV). Default: %(default)s",
    )
    p.add_argument(
        "--outdir",
        default="pathogenic_variants",
        help="Base output dir for per-case VCFs. Default: %(default)s",
    )
    p.add_argument(
        "--sleep",
        type=float,
        default=0.2,
        help="Sleep seconds between REST calls to be polite. Default: %(default)s",
    )
    p.add_argument(
        "--dedup-info",
        action="store_true",
        help="When duplicate coordinates are found, append all HGVS notations "
             "to the INFO field (HGVS=first;HGVS_ALT2=second;...). "
             "Without this flag, only the first HGVS is kept.",
    )
    return p.parse_args()


def load_hgvs_by_case(merged_path):
    """
    读取 merged 表，按 (Case, Reference_genome) 分组 HGVS 列表。

    列结构：
      #Case   Gene   cHGVS   pHGVS   Ref seq   Phenotype   Disease   Reference_genome
    """
    if not os.path.isfile(merged_path):
        sys.stderr.write(f"[ERROR] merged file not found: {merged_path}\n")
        sys.exit(1)

    by_case = defaultdict(list)

    with open(merged_path, newline="") as f:
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

            # 一个 row 里可能有多个 cHGVS，用逗号分开
            for c in chgvs_raw.split(","):
                c = c.strip()
                if not c or c == "-":
                    continue
                hgvs = f"{refseq}:{c}"
                by_case[(case, genome)].append(hgvs)

    # 去重（按 HGVS 字符串）
    for key in list(by_case.keys()):
        uniq = sorted(set(by_case[key]))
        by_case[key] = uniq

    return by_case


def pick_server_for_genome(genome):
    """
    根据 Reference_genome 选择 Ensembl REST server。
    """
    if genome == "GRCh37":
        return SERVER_37
    elif genome == "GRCh38":
        return SERVER_38
    else:
        return None


def chunked(iterable, size):
    """简单的分块函数。"""
    it = iter(iterable)
    while True:
        chunk = []
        try:
            for _ in range(size):
                chunk.append(next(it))
        except StopIteration:
            if chunk:
                yield chunk
            break
        yield chunk


def query_ensembl_rest(server, hgvs_list):
    """
    调用 Ensembl REST VEP HGVS endpoint，返回 JSON。
    """
    url = server + EXT_HGVS
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
    }
    payload = {"hgvs_notations": hgvs_list}

    r = requests.post(url, headers=headers, data=json.dumps(payload))
    if not r.ok:
        sys.stderr.write(
            f"[ERROR] REST call failed ({r.status_code}) {url}\n{r.text}\n"
        )
        r.raise_for_status()
    return r.json()


def dedup_by_coordinate(records, keep_all_hgvs=False):
    """
    按 (CHROM, POS, REF, ALT) 去重。

    Args:
        records: list of dicts with keys: chrom, pos, ref, alt, hgvs
        keep_all_hgvs: if True, 当遇到重复坐标时，将后续 HGVS 追加到
                       记录的 extra_hgvs 列表中

    Returns:
        deduped list of dicts; each dict may have an extra 'extra_hgvs' key
    """
    seen = OrderedDict()  # key → record (preserves insertion order)
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
        print(f"[DEDUP] Removed {dup_count} duplicate coordinate(s)")

    return list(seen.values())


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def write_case_vcf(case, genome, records, outdir):
    """
    写出一个 case 的 VCF 文件。

    records: list of dicts with keys: chrom, pos, ref, alt, hgvs
             optionally: extra_hgvs (list of additional HGVS notations)
    """
    subdir = os.path.join(outdir, genome)
    ensure_dir(subdir)
    out_path = os.path.join(subdir, f"case{case}.vcf")

    with open(out_path, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##source=EnsemblREST_VEP_HGVS_dedup\n")
        out.write(f"##reference_genome={genome}\n")
        out.write("##INFO=<ID=HGVS,Number=1,Type=String,Description=\"Original HGVS notation\">\n")
        out.write("##INFO=<ID=HGVS_ALT,Number=.,Type=String,Description=\"Additional HGVS notations at same coordinate (deduped)\">\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for rec in records:
            chrom = rec["chrom"]
            pos = rec["pos"]
            ref = rec["ref"]
            alt = rec["alt"]
            hgvs = rec["hgvs"].replace(" ", "_")  # 避免 INFO 里有空格
            info_fields = [f"HGVS={hgvs}"]

            # 如果有额外 HGVS（来自同一坐标的重复记录）
            if "extra_hgvs" in rec and rec["extra_hgvs"]:
                alt_hgvs = ",".join(h.replace(" ", "_") for h in rec["extra_hgvs"])
                info_fields.append(f"HGVS_ALT={alt_hgvs}")

            info = ";".join(info_fields)
            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{info}\n")

    print(f"[OK] Case {case} ({genome}), variants: {len(records)}, VCF: {out_path}")


def main():
    args = parse_args()
    merged_path = os.path.abspath(args.merged)
    outdir = os.path.abspath(args.outdir)

    print(f"[INFO] Merged table : {merged_path}")
    print(f"[INFO] Output dir   : {outdir}")
    print(f"[INFO] Dedup info   : {args.dedup_info}")
    print()

    by_case = load_hgvs_by_case(merged_path)
    print(f"[INFO] Total (Case,Genome) groups: {len(by_case)}")

    total_dup = 0

    for (case, genome), hgvs_list in sorted(by_case.items(), key=lambda x: int(x[0][0])):
        server = pick_server_for_genome(genome)
        if server is None:
            sys.stderr.write(
                f"[WARN] Case {case}: unsupported Reference_genome '{genome}', skip.\n"
            )
            continue

        print(f"\n=== Case {case} / {genome} ===")
        print(f"[INFO] HGVS count: {len(hgvs_list)}")
        all_records = []

        for chunk in chunked(hgvs_list, MAX_PER_POST):
            print(f"[INFO]  Querying {len(chunk)} HGVS to {server} ...")
            try:
                data = query_ensembl_rest(server, chunk)
            except Exception as e:
                sys.stderr.write(
                    f"[ERROR] REST query failed for Case {case} ({genome}): {e}\n"
                )
                continue

            # data 是一个数组，每个元素对应一个 HGVS 的结果
            for item in data:
                v_input = item.get("input")
                vcf_string = item.get("vcf_string")  # 形如 "1-1196415-T-C"
                if not vcf_string:
                    sys.stderr.write(
                        f"[WARN] Missing vcf_string for input '{v_input}', skip.\n"
                    )
                    continue

                parts = vcf_string.split("-")
                if len(parts) != 4:
                    sys.stderr.write(
                        f"[WARN] Unexpected vcf_string format '{vcf_string}' for input '{v_input}', skip.\n"
                    )
                    continue

                chrom, pos, ref, alt = parts
                rec = {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "hgvs": v_input,
                }
                all_records.append(rec)

            time.sleep(args.sleep)

        if not all_records:
            sys.stderr.write(
                f"[WARN] Case {case} ({genome}): no valid records from Ensembl, skip writing VCF.\n"
            )
            continue

        # ===== 新增：按坐标去重 =====
        before_dedup = len(all_records)
        all_records = dedup_by_coordinate(all_records, keep_all_hgvs=args.dedup_info)
        after_dedup = len(all_records)
        dup_this_case = before_dedup - after_dedup
        total_dup += dup_this_case
        if dup_this_case > 0:
            print(f"[DEDUP] Case {case} ({genome}): {before_dedup} → {after_dedup} records "
                  f"({dup_this_case} duplicates removed by coordinate)")
        # ===== 去重结束 =====

        write_case_vcf(case, genome, all_records, outdir)

    print(f"\n[INFO] All done. Total coordinate duplicates removed: {total_dup}")


if __name__ == "__main__":
    main()
