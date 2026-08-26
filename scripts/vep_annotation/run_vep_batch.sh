#!/bin/bash
# VEP批量标注脚本 (GRCh37) - 并行版本
# 用法: ./run_vep_batch.sh [并行数]
# 示例: ./run_vep_batch.sh 4  (并行运行4个VEP任务)

set -e

# 配置
VCF_DIR="genotype/GRCh37"
VEP_DATA_DIR="$(pwd)/vep_data"
THREADS=8                              # 每个VEP任务的线程数
PARALLEL_JOBS=${1:-4}                  # 并行任务数，默认4
LOG_DIR="$(pwd)/vep_logs"

# 计算总资源使用
TOTAL_THREADS=$((THREADS * PARALLEL_JOBS))
echo "========================================"
echo "VEP批量标注 (GRCh37) - 并行模式"
echo "========================================"
echo "并行任务数: $PARALLEL_JOBS"
echo "每任务线程: $THREADS"
echo "总线程占用: $TOTAL_THREADS / $(nproc)"
echo "========================================"

# 检查输入目录
if [ ! -d "$VCF_DIR" ]; then
    echo "错误: VCF目录不存在: $VCF_DIR"
    exit 1
fi

# 创建日志目录
mkdir -p "$LOG_DIR"

# 统计信息
total_count=$(find "$VCF_DIR" -maxdepth 1 -name "*.vcf.gz" ! -name "*.vep.vcf.gz" ! -name "*.tbi" | wc -l)
completed_count=$(find "$VCF_DIR" -maxdepth 1 -name "*.vep.vcf.gz" | wc -l)
echo "总VCF文件: $total_count"
echo "已完成: $completed_count"
echo "待处理: $((total_count - completed_count))"
echo "日志目录: $LOG_DIR"
echo "========================================"

# 并行控制
declare -a pids=()
declare -A pid_to_file=()
running=0
processed=0
skipped=0
failed=0
queue=()

# 单个文件处理函数
process_vcf() {
    local vcf="$1"
    local filename=$(basename "$vcf")
    local vep_output="${vcf%.vcf.gz}.vep.vcf.gz"
    local log_file="$LOG_DIR/${filename%.vcf.gz}.log"

    # 检查输出文件是否已存在且完整
    if [ -f "$vep_output" ]; then
        local file_size=$(stat -c%s "$vep_output" 2>/dev/null || echo "0")
        if [ "$file_size" -gt 100000 ]; then
            echo "[跳过] $filename" >> "$LOG_DIR/summary.log"
            return 0  # skipped
        else
            rm -f "$vep_output"
        fi
    fi

    local docker_input="/data/$filename"
    local docker_output="/data/${filename%.vcf.gz}.vep.vcf.gz"

    echo "[$(date '+%H:%M:%S')] 开始: $filename" | tee -a "$log_file"

    if docker run --rm \
        -v "$VEP_DATA_DIR:/vep_data" \
        -v "$(pwd)/$VCF_DIR:/data" \
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
        --plugin MaxEntScan,/vep_data/fordownload,SWA,NCSS \
        >> "$log_file" 2>&1; then

        if [ -f "$vep_output" ] && [ -s "$vep_output" ]; then
            local out_size=$(stat -c%s "$vep_output")
            echo "[$(date '+%H:%M:%S')] ✓ 完成: $filename (${out_size} bytes)" | tee -a "$log_file"
            return 1  # success
        else
            echo "[$(date '+%H:%M:%S')] ✗ 失败: $filename (输出为空)" | tee -a "$log_file"
            rm -f "$vep_output"
            return 2  # failed
        fi
    else
        echo "[$(date '+%H:%M:%S')] ✗ 失败: $filename (VEP错误)" | tee -a "$log_file"
        rm -f "$vep_output"
        return 2  # failed
    fi
}

# 收集完成的任务
collect_finished() {
    local new_pids=()
    for pid in "${pids[@]}"; do
        if ! kill -0 "$pid" 2>/dev/null; then
            wait "$pid"
            local result=$?
            case $result in
                0) ((skipped++)) ;;
                1) ((processed++)) ;;
                2) ((failed++)) ;;
            esac
            ((running--))
        else
            new_pids+=("$pid")
        fi
    done
    pids=("${new_pids[@]}")
}

echo "[$(date '+%H:%M:%S')] 开始处理..."

# 主循环
for vcf in "$VCF_DIR"/*.vcf.gz; do
    [ -e "$vcf" ] || continue
    filename=$(basename "$vcf")
    [[ "$filename" == *.vep.vcf.gz ]] && continue

    # 等待有空闲槽位
    while [ $running -ge $PARALLEL_JOBS ]; do
        sleep 5
        collect_finished
    done

    # 启动新任务
    process_vcf "$vcf" &
    pids+=($!)
    pid_to_file[$!]="$filename"
    ((running++))

    echo -ne "\r运行中: $running / $PARALLEL_JOBS | 已完成: $processed | 跳过: $skipped | 失败: $failed"
done

# 等待所有任务完成
echo ""
echo "[$(date '+%H:%M:%S')] 等待剩余任务完成..."
while [ $running -gt 0 ]; do
    collect_finished
    echo -ne "\r运行中: $running | 已完成: $processed | 跳过: $skipped | 失败: $failed"
    sleep 5
done

# 最终统计
echo ""
echo ""
echo "========================================"
echo "批量标注完成! $(date '+%Y-%m-%d %H:%M:%S')"
echo "========================================"
echo "本次处理: $processed"
echo "跳过(已存在): $skipped"
echo "失败: $failed"
echo "日志目录: $LOG_DIR"
echo "========================================"
