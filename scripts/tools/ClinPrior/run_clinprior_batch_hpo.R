# 文件名: run_clinprior_batch_hpo.R
# 作用: 自动修复数据库并批量运行 ClinPrior HPO-only 基因排序 (GRCh37 + GRCh38)
#
# 生成 clinprior_results_hpo/ 目录:
#   每个 HPO 文件 case{N}_{SampleID}_hpos.txt -> case{N}_{SampleID}_hpos_GRCh{37|38}_result.csv
#   内容为 proteinScore() + MatrixPropagation(alpha=0.2) 的原始输出,
#   按 PriorFunct 降序排列, 列: Symbol, geneID, PriorFunct, Symbol, geneID, PriorPhys
#
# 说明: proteinScore() 要求 >=2 个 HPO term; 单 term 案例会触发 tryCatch 被跳过 (不写文件)。

library(ClinPrior)
library(piggyback) # 用于自动下载缺失数据

# ==============================================================================
# 1. 基础配置
# ==============================================================================
# 运行前设置 GENETOOLS_BASE 环境变量指向项目根目录, 例如:
#   Sys.setenv(GENETOOLS_BASE = "/path/to/genetools")
# 或在 shell 中:  export GENETOOLS_BASE=/path/to/genetools
base_dir <- Sys.getenv("GENETOOLS_BASE", unset = "")
if (nchar(base_dir) == 0) {
  stop("Set the GENETOOLS_BASE environment variable to your project root ",
       "(e.g. Sys.setenv(GENETOOLS_BASE='/path/to/genetools') or ",
       "export GENETOOLS_BASE=/path/to/genetools).")
}

# 输入文件夹
dir_37 <- file.path(base_dir, "phenotype/GRCh37")
dir_38 <- file.path(base_dir, "phenotype/GRCh38")

# 统一输出文件夹 (HPO-only 基因排序结果)
output_dir <- file.path(base_dir, "clinprior_results_hpo")

# 数据库根目录 (通常在 /usr/local/lib/R/site-library/ClinPrior/extdata)
db_root <- system.file("extdata", package = "ClinPrior")

# 自动创建输出目录
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("==========================================================\n")
cat(">>> ClinPrior 全自动批量分析脚本启动\n")
cat(">>> 输出目录:", output_dir, "\n")
cat("==========================================================\n")


# ==============================================================================
# 2. 定义辅助函数: 自动检查并修复数据库
# ==============================================================================
check_and_repair_db <- function(folder_name, tag_name) {
  target_path <- file.path(db_root, folder_name)

  # 1. 如果文件夹不存在，创建它
  if (!dir.exists(target_path)) dir.create(target_path, recursive = TRUE)

  # 2. 检查核心文件是否存在
  files <- list.files(target_path)
  has_core_files <- any(grepl("Physical_interactome.RData", files)) &&
                    any(grepl("Functional_interactome.RData", files))

  if (!has_core_files) {
    cat(paste0("\n⚠️  检测到 ", folder_name, " 缺少核心文件，正在自动下载修复 (这可能需要几分钟)...\n"))
    tryCatch({
      pb_download(repo = "aschluter/ClinPrior",
                  tag = tag_name,
                  dest = target_path,
                  overwrite = TRUE)
      cat("✅ 下载完成！\n")
    }, error = function(e) {
      cat("❌ 下载失败，请检查网络连接。\n")
      stop(e)
    })
  } else {
    cat(paste0("✅ ", folder_name, " 数据库完整。\n"))
  }
}


# ==============================================================================
# 3. 定义核心处理函数
# ==============================================================================
process_batch <- function(input_path, version_name, db_folder) {

  cat(paste0("\n----------------------------------------------------------\n"))
  cat(paste0(">>> [阶段开始] 正在处理 ", version_name, " 数据...\n"))

  # A. 加载数据库
  db_path <- file.path(db_root, db_folder)
  cat(paste0(">>> 正在加载数据库: ", db_path, "\n"))

  tryCatch({
    # 必须加载到全局环境，否则 MatrixPropagation 找不到
    load(file.path(db_path, "Physical_interactome.RData"), envir = .GlobalEnv)
    load(file.path(db_path, "Functional_interactome.RData"), envir = .GlobalEnv)
  }, error = function(e) {
    cat("❌ 严重错误: 无法加载 .RData 文件。请确保上一步修复成功。\n")
    return(NULL)
  })

  # B. 获取文件列表
  files <- list.files(input_path, pattern = "\\.txt$", full.names = TRUE)
  if (length(files) == 0) {
    cat(paste0("⚠️  警告: ", input_path, " 下没有 .txt 文件，跳过此批次。\n"))
    return(NULL)
  }

  # C. 循环处理
  count <- 0
  total <- length(files)

  for (f in files) {
    fname <- basename(f)
    # 简单的进度提示
    cat(sprintf("   [%d/%d] 处理: %s ... ", count + 1, total, fname))

    tryCatch({
      # 1. 读取数据
      dat <- read.table(f, header = FALSE, stringsAsFactors = FALSE)
      hpos <- unique(dat[, 1])

      # 2. 核心计算
      Y <- proteinScore(hpos)
      res <- MatrixPropagation(Y, alpha=0.2)

      # 3. 格式转换与排序 (关键步骤)
      df_res <- as.data.frame(res)

      # 优先按 PriorFunct (功能分) 降序排列
      if ("PriorFunct" %in% colnames(df_res)) {
        df_sorted <- df_res[order(df_res$PriorFunct, decreasing = TRUE), ]
      } else {
        # 如果列名不匹配，按第一列排
        df_sorted <- df_res[order(df_res[, 1], decreasing = TRUE), ]
      }

      # 4. 保存文件 (添加版本后缀)
      out_name <- paste0(tools::file_path_sans_ext(fname), "_", version_name, "_result.csv")
      out_file <- file.path(output_dir, out_name)

      write.csv(df_sorted, out_file, row.names = FALSE)
      cat("完成\n")
      count <- count + 1

    }, error = function(e) {
      cat(paste0("失败! (", e$message, ")\n"))
    })
  }
  cat(paste0(">>> ", version_name, " 处理完毕。\n"))
}


# ==============================================================================
# 4. 正式执行 (Step-by-Step)
# ==============================================================================

# --- 第一步: 检查并修复数据库 (解决缺文件问题) ---
cat("\n>>> Step 1: 检查数据库完整性...\n")
check_and_repair_db("assembly37", "v2.0-assemblyGRCh37")
check_and_repair_db("assembly38", "v2.0-assemblyGRCh38")

# --- 第二步: 运行 GRCh37 批次 ---
# 注意: 这里明确指定使用 "assembly37" 文件夹
if (dir.exists(dir_37)) {
  process_batch(input_path = dir_37,
                version_name = "GRCh37",
                db_folder = "assembly37")
} else {
  cat(paste0("\n⚠️  输入目录不存在: ", dir_37, "\n"))
}

# --- 第三步: 运行 GRCh38 批次 ---
if (dir.exists(dir_38)) {
  process_batch(input_path = dir_38,
                version_name = "GRCh38",
                db_folder = "assembly38")
} else {
  cat(paste0("\n⚠️  输入目录不存在: ", dir_38, "\n"))
}

cat("\n==========================================================\n")
cat("🎉 所有任务已结束！结果保存在:", output_dir, "\n")
