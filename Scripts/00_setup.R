# ============================================================
# 00_setup.R
# 安装并加载所有必需的R包，并设置固定的工作目录
# ============================================================

# ---- 0. 设置项目工作目录 (请根据实际情况修改) ----
# 使用绝对路径确保无论在何处启动R，都能定位到项目文件夹
project_dir <- "/Users/bing/DH"
if (!dir.exists(project_dir)) {
  dir.create(project_dir, recursive = TRUE)
}
setwd(project_dir)
cat("当前工作目录已设置为:", getwd(), "\n")

# 创建一个专用输出文件夹，用于存放所有结果
output_dir <- file.path(project_dir, "output")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("已创建输出文件夹:", output_dir, "\n")
}

# ---- 1. 安装并加载所有必需的R包 ----
required_packages <- c(
  "nhanesA",      # NHANES数据下载
  "survival",     # 生存分析
  "rms",          # Nomogram和cph模型
  "timeROC",      # 时间依赖ROC
  "survminer",    # KM曲线可视化
  "ggplot2",      # 绘图
  "dplyr",        # 数据操作
  "broom",        # 模型整洁输出
  "nricens",      # 净重分类改善
  "haven",        # 读取XPT文件（备选）
  "TwoSampleMR",  # 孟德尔随机化
  "MRPRESSO",     # MR多效性检测
  "ieugwasr"      # GWAS数据接口（依赖）
)

# 安装缺失的包
installed <- rownames(installed.packages())
for(pkg in required_packages){
  if(!pkg %in% installed){
    install.packages(pkg, dependencies = TRUE)
  }
}

# 加载所有包
lapply(required_packages, library, character.only = TRUE)

# 如果ieugwasr需要单独安装，使用remotes
if(!require(ieugwasr)){
  install.packages("remotes")
  remotes::install_github("MRCIEU/ieugwasr")
}

message("All packages installed and loaded successfully.")
