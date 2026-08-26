#!/bin/bash
# VEP单文件标注脚本 (用于Condor)
# 用法: ./run_vep_single.sh <vcf_file>

set -e

if [ -z "$1" ]; then
    echo "用法: $0 <vcf_filename>"
    exit 1
fi

VCF_FILE="$1"
# Base directory of the benchmark + annotation data.
# Set GENETOOLS_BASE before running, e.g.:
#   export GENETOOLS_BASE=/path/to/genetools
if [ -z "${GENETOOLS_BASE:-}" ]; then
    echo "错误: 请先设置 GENETOOLS_BASE 环境变量 (项目根目录)"
    echo "  例如: export GENETOOLS_BASE=/path/to/genetools"
    exit 1
fi
BASE_DIR="$GENETOOLS_BASE"
VCF_DIR="$BASE_DIR/genotype/GRCh37"
VEP_DATA_DIR="$BASE_DIR/vep_data"
THREADS=4

# 检查docker是否可用
if ! command -v docker &> /dev/null; then
    # 尝试常见的docker路径
    for docker_path in /usr/bin/docker /usr/local/bin/docker /opt/docker/bin/docker; do
        if [ -x "$docker_path" ]; then
            export PATH="$PATH:$(dirname $docker_path)"
            break
        fi
    done
    if ! command -v docker &> /dev/null; then
        echo "错误: docker未安装或不在PATH中"
        exit 1
    fi
fi

# 转换为绝对路径
VCF_PATH="$VCF_DIR/$VCF_FILE"
VEP_OUTPUT="${VCF_PATH%.vcf.gz}.vep.vcf.gz"

# 检查输入文件
if [ ! -f "$VCF_PATH" ]; then
    echo "错误: 输入文件不存在: $VCF_PATH"
    exit 1
fi

# 检查输出是否已存在且完整
if [ -f "$VEP_OUTPUT" ]; then
    file_size=$(stat -c%s "$VEP_OUTPUT" 2>/dev/null || echo "0")
    if [ "$file_size" -gt 100000 ]; then
        echo "跳过: $VCF_FILE (已有完整输出)"
        exit 0
    else
        echo "重试: $VCF_FILE (输出不完整)"
        rm -f "$VEP_OUTPUT"
    fi
fi

echo "开始处理: $VCF_FILE"
echo "主机: $(hostname)"
echo "时间: $(date)"

docker_input="/data/$VCF_FILE"
docker_output="/data/${VCF_FILE%.vcf.gz}.vep.vcf.gz"

docker run --rm \
    -v "$VEP_DATA_DIR:/vep_data" \
    -v "$VCF_DIR:/data" \
    ensemblorg/ensembl-vep:latest \
    vep \
    --dir_cache /vep_data \
    --dir_plugins /opt/vep/.vep/Plugins \
    --assembly GRCh37 \
    --fasta /vep_data/homo_sapiens/115_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
    --input_file "$docker_input" \
    --output_file "$docker_output" \
    --compress_output gzip \
    --af --af_gnomade --af_gnomadg --appris --biotype \
    --buffer_size 10000 --cache --ccds --check_existing \
    --distance 5000 --filter_common --force --fork $THREADS \
    --hgvs --mane --no_stats --offline \
    --pick_allele --polyphen b --quiet \
    --regulatory --safe --sift b \
    --symbol --transcript_version --tsl --var_synonyms --vcf \
    --plugin CADD,snv=/vep_data/CADD/GRCh37/whole_genome_SNVs.tsv.gz \
    --plugin MaxEntScan,/vep_data/fordownload,SWA,NCSS

# 验证输出
if [ -f "$VEP_OUTPUT" ] && [ -s "$VEP_OUTPUT" ]; then
    out_size=$(stat -c%s "$VEP_OUTPUT")
    echo "完成: $VCF_FILE (${out_size} bytes)"
    echo "时间: $(date)"
    exit 0
else
    echo "失败: $VCF_FILE (输出为空)"
    rm -f "$VEP_OUTPUT"
    exit 1
fi
