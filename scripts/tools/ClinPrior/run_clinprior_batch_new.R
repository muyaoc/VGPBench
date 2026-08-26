# ==============================================================================
# ClinPrior 批量处理脚本 (New - 遵循 ClinPrior README API)
# 使用 readVCF() + priorBestVariantVcfR() 标准流程
#
# 与旧脚本的区别:
#   1. 使用 ClinPrior 标准 API (readVCF/priorBestVariantVcfR)
#   2. GT 转换仅替换样本列中的 |, 保留 CSQ 字段
#   3. 遵循 README 参数: geneQuality=20, readDepth=10,
#      distSplicThreshold="full"(→Inf), synonymous=TRUE,
#      discard_zero_score=TRUE
#   4. 输出: clinprior_results_new_<date>/<case#>_<sample_id>_VariantPrior.csv
# ==============================================================================

library(ClinPrior)
library(vcfR)

# ==============================================================================
# 配置路径
# ==============================================================================
# Set GENETOOLS_BASE before running, e.g.:
#   Sys.setenv(GENETOOLS_BASE = "/path/to/genetools")
# or from the shell:  export GENETOOLS_BASE=/path/to/genetools
base_dir <- Sys.getenv("GENETOOLS_BASE", unset = "")
if (nchar(base_dir) == 0) {
  stop("Set the GENETOOLS_BASE environment variable to your project root ",
       "(e.g. Sys.setenv(GENETOOLS_BASE='/path/to/genetools') or ",
       "export GENETOOLS_BASE=/path/to/genetools).")
}

vcf_dir <- file.path(base_dir, "genotype/GRCh37")
hpo_dir <- file.path(base_dir, "phenotype/GRCh37")

date_str <- format(Sys.Date(), "%Y%m%d")
output_dir <- file.path(base_dir, paste0("clinprior_results_new_", date_str))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

temp_dir <- file.path(output_dir, ".temp_unphased")
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(output_dir, "batch_log.txt")

# ==============================================================================
# 命令行参数 (可选)
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
skip_existing <- TRUE
case_subset <- NULL

if (length(args) > 0) {
  for (arg in args) {
    if (arg == "--force") {
      skip_existing <- FALSE
    } else if (arg == "--help") {
      cat("Usage: Rscript run_clinprior_batch_new.R [--force] [case1 case2 ...]\n")
      cat("  --force   Overwrite existing results\n")
      cat("  case#     Process only specified cases (e.g., case1_HG00407)\n")
      quit(status = 0)
    } else {
      case_subset <- c(case_subset, arg)
    }
  }
}

# ==============================================================================
# GT 转换: 仅替换样本列(col>=10)中的 |, 保留 INFO/CSQ 字段
# ==============================================================================
convertPhasedGT <- function(input_vcf, output_vcf) {
  cmd <- paste0(
    "zcat '", input_vcf, "' | ",
    "awk -F'\\t' -v OFS='\\t' '{for(i=10;i<=NF;i++) gsub(/\\|/,\"/\",$i); print}' | ",
    "gzip > '", output_vcf, "'"
  )
  result <- system(cmd, intern = FALSE)
  if (result != 0) {
    stop(paste("GT conversion failed for", input_vcf))
  }
  return(output_vcf)
}

# ==============================================================================
# 单样本处理函数 - 遵循 ClinPrior README API
# ==============================================================================
processSample <- function(case_info, sample_name, vcf_file, hpo_file) {
  output_csv <- file.path(output_dir, paste0(case_info, "_VariantPrior.csv"))

  if (skip_existing && file.exists(output_csv) && file.size(output_csv) > 100) {
    cat(paste0("[SKIP] ", case_info, " already done\n"))
    return(TRUE)
  }

  cat(paste0("\n=== ", case_info, " ===\n"))

  tryCatch({
    # ------------------------------------------------------------------
    # Step 1: HPO -> GlobalPhenotypicScore
    # README: Y <- proteinScore(HPOpatient)
    #         ClinPriorGeneScore <- MatrixPropagation(Y, alpha = 0.2)
    # ------------------------------------------------------------------
    hpo_list <- unique(read.csv(hpo_file, header = FALSE, sep = "\t")[[1]])
    hpo_list <- as.character(trimws(hpo_list))
    hpo_list <- hpo_list[hpo_list != ""]

    if (length(hpo_list) < 2) {
      cat(paste0("  Only ", length(hpo_list), " HPO term(s), duplicating to meet minimum\n"))
      hpo_list <- rep(hpo_list, 3)
    }

    cat(paste0("  HPO terms: ", length(unique(hpo_list)), "\n"))

    Y <- proteinScore(hpo_list)
    GlobalPhenotypicScore <- MatrixPropagation(Y, alpha = 0.2)

    # ------------------------------------------------------------------
    # Step 2: 转换 phased GT (0|1 -> 0/1)
    # ClinPrior readVCF 使用 grep("0/0",...) 过滤, 只识别 unphased 格式
    # ------------------------------------------------------------------
    unphased_vcf <- file.path(temp_dir, paste0(case_info, ".unphased.vcf.gz"))
    if (!file.exists(unphased_vcf)) {
      cat("  Converting phased GT -> unphased...\n")
      convertPhasedGT(vcf_file, unphased_vcf)
    }

    # ------------------------------------------------------------------
    # Step 3: 读取 VCF
    # README: variants <- readVCF(sampleName, vcfFile, assembly, ...)
    # 安装版本: readVCF(sampleName, variants=vcfR_obj, assembly, ...)
    # ------------------------------------------------------------------
    cat("  Reading VCF...\n")
    variants <- read.vcfR(unphased_vcf, verbose = FALSE)
    cat(paste0("  Total variants: ", nrow(variants), "\n"))

    # ------------------------------------------------------------------
    # Step 4: 过滤变异 (readVCF - ClinPrior README API)
    # README参数:
    #   geneQuality        = 20
    #   readDepth          = 10
    #   distSplicThreshold = "full" (→ Inf, 保留所有内含子变异)
    #   synonymous         = TRUE
    # ------------------------------------------------------------------
    cat("  Filtering with readVCF...\n")
    variantsFiltered <- readVCF(
      sampleName         = sample_name,
      variants           = variants,
      assembly           = "assembly37",
      geneQuality        = 20,
      readDepth          = 10,
      distSplicThreshold = Inf,
      synonymous         = TRUE
    )
    cat(paste0("  After filter: ", nrow(variantsFiltered), " variants\n"))

    if (nrow(variantsFiltered) == 0) {
      cat("  No variants after filtering\n")
      write.csv(data.frame(), file = output_csv, row.names = FALSE)
      return(TRUE)
    }

    # ------------------------------------------------------------------
    # Step 5: 变异优先级排序 (priorBestVariantVcfR - ClinPrior README API)
    # README: results <- priorBestVariant(variants, sampleName,
    #             GlobalPhenotypicScore, assembly, splicingMode="vep",
    #             discard_zero_score=TRUE)
    # 安装版本: priorBestVariantVcfR (splicingMode已默认使用VEP MaxEntScan)
    # ------------------------------------------------------------------
    cat("  Prioritizing variants...\n")
    results <- priorBestVariantVcfR(
      variants              = variantsFiltered,
      sampleName            = sample_name,
      GlobalPhenotypicScore = GlobalPhenotypicScore,
      assembly              = "assembly37"
    )

    # discard_zero_score = TRUE (README)
    if ("ClinPriorScore" %in% colnames(results)) {
      n_before <- nrow(results)
      results <- results[results$ClinPriorScore > 0, ]
      cat(paste0("  Discarded ", n_before - nrow(results), " zero-score variants\n"))
    }

    # 重新编号 ClinPriorPosition
    if (nrow(results) > 0 && "ClinPriorPosition" %in% colnames(results)) {
      results$ClinPriorPosition <- 1:nrow(results)
    }

    cat(paste0("  Results: ", nrow(results), " candidate variants\n"))

    if (nrow(results) > 0) {
      cat("  Top variant:\n")
      print(head(results[, intersect(
        c("ClinPriorPosition", "CHROM", "POS", "REF", "ALT",
          "genesList", "clinvar", "knownDisease", "Consequence",
          "cDNA", "Protein", "ClinPriorScore"),
        colnames(results)
      ), drop = FALSE], 1))
    }

    # ------------------------------------------------------------------
    # Step 6: 保存结果
    # 仅保留 README 核心列: ClinPriorPosition, CHROM, POS, REF, ALT,
    #   genesList, clinvar, knownDisease, Consequence, cDNA, Protein,
    #   ClinPriorScore
    # ------------------------------------------------------------------
    readme_cols <- c("ClinPriorPosition", "CHROM", "POS", "REF", "ALT",
                     "genesList", "clinvar", "knownDisease", "Consequence",
                     "cDNA", "Protein", "ClinPriorScore")
    cols_to_keep <- intersect(readme_cols, colnames(results))
    results <- results[, cols_to_keep, drop = FALSE]

    write.csv(results, file = output_csv, row.names = FALSE)
    cat(paste0("  Saved: ", output_csv, "\n"))

    # 释放内存
    rm(variants, variantsFiltered, results, GlobalPhenotypicScore)
    gc()

    return(TRUE)

  }, error = function(e) {
    cat(paste0("  ERROR: ", e$message, "\n"))
    # 写入日志
    cat(paste0(Sys.time(), " | ", case_info, " | ERROR | ", e$message, "\n"),
        file = log_file, append = TRUE)
    return(FALSE)
  })
}

# ==============================================================================
# 主程序
# ==============================================================================
cat("============================================================\n")
cat("ClinPrior Batch Analysis (New - README API)\n")
cat(paste0("Date: ", date_str, "\n"))
cat(paste0("Output: ", output_dir, "\n"))
cat(paste0("Skip existing: ", skip_existing, "\n"))
cat("============================================================\n\n")

# 写入日志头
cat(paste0(Sys.time(), " | START | VCF_dir=", vcf_dir, "\n"),
    file = log_file, append = FALSE)

# 构建文件映射
vcf_files <- list.files(vcf_dir, pattern = "\\.vep\\.vcf\\.gz$", full.names = TRUE)
hpo_files <- list.files(hpo_dir, pattern = "_hpos\\.txt$", full.names = TRUE)

cat(paste0("VEP VCF files: ", length(vcf_files), "\n"))
cat(paste0("HPO files: ", length(hpo_files), "\n"))

# case_info -> vcf_file
vcf_map <- list()
for (f in vcf_files) {
  m <- regmatches(basename(f), regexpr("case[0-9]+_[A-Za-z0-9]+", basename(f)))
  if (length(m) > 0) vcf_map[[m]] <- f
}

# case_info -> hpo_file
hpo_map <- list()
for (f in hpo_files) {
  m <- regmatches(basename(f), regexpr("case[0-9]+_[A-Za-z0-9]+", basename(f)))
  if (length(m) > 0) hpo_map[[m]] <- f
}

# 需要同时有 VCF 和 HPO 的 case
common <- intersect(names(vcf_map), names(hpo_map))

# 仅 VCF 没有 HPO 的 case
vcf_only <- setdiff(names(vcf_map), names(hpo_map))

cat(paste0("Cases with VCF+HPO: ", length(common), "\n"))
cat(paste0("Cases with VCF only (no HPO, skipped): ", length(vcf_only), "\n\n"))

if (length(vcf_only) > 0 && length(vcf_only) <= 10) {
  cat(paste0("  Missing HPO: ", paste(vcf_only, collapse = ", "), "\n\n"))
}

# 如果指定了子集, 只处理指定的 case
if (!is.null(case_subset)) {
  common <- intersect(common, case_subset)
  cat(paste0("Processing subset: ", length(common), " cases\n\n"))
}

# 排序处理
common <- sort(common)

# ==============================================================================
# 逐样本处理
# ==============================================================================
n_ok <- 0
n_fail <- 0
n_skip <- 0

t_start <- Sys.time()

for (i in seq_along(common)) {
  ci <- common[i]
  sn <- sub("case[0-9]+_", "", ci)

  # 进度显示
  if (i %% 50 == 0 || i == length(common)) {
    elapsed <- difftime(Sys.time(), t_start, units = "mins")
    rate <- i / as.numeric(elapsed)
    eta <- (length(common) - i) / rate
    cat(paste0("\n--- Progress: ", i, "/", length(common),
               " (", round(rate, 1), " cases/min, ETA ~", round(eta, 0), " min) ---\n"))
  }

  ok <- processSample(ci, sn, vcf_map[[ci]], hpo_map[[ci]])

  if (ok) {
    if (file.exists(file.path(output_dir, paste0(ci, "_VariantPrior.csv")))) {
      n_ok <- n_ok + 1
    } else {
      n_skip <- n_skip + 1
    }
  } else {
    n_fail <- n_fail + 1
  }
}

# ==============================================================================
# 汇总
# ==============================================================================
elapsed_total <- round(difftime(Sys.time(), t_start, units = "mins"), 1)

cat("\n============================================================\n")
cat("Summary\n")
cat("============================================================\n")
cat(paste0("Total:     ", length(common), " cases\n"))
cat(paste0("Success:   ", n_ok, "\n"))
cat(paste0("Failed:    ", n_fail, "\n"))
cat(paste0("Skipped:   ", n_skip, "\n"))
cat(paste0("Elapsed:   ", elapsed_total, " min\n"))
cat(paste0("Output:    ", output_dir, "\n"))
cat(paste0("Log:       ", log_file, "\n"))

# 列出失败样本
if (n_fail > 0) {
  cat("\nFailed cases (see log for details):\n")
  if (file.exists(log_file)) {
    failed_lines <- grep("ERROR", readLines(log_file), value = TRUE)
    for (l in failed_lines) cat(paste0("  ", l, "\n"))
  }
}

cat("============================================================\n")

# 写入日志
cat(paste0(Sys.time(), " | DONE | success=", n_ok, " failed=", n_fail,
           " skipped=", n_skip, " elapsed=", elapsed_total, "min\n"),
    file = log_file, append = TRUE)
