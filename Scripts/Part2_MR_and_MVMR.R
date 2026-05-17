# ==============================================================================
# 项目名称: Double Heart (DH)
# 分析模块: 第二部分：双向与多因素孟德尔随机化 (MR) 分析
# 步骤名称: Step 1: 心理躯体症状与心血管结局的高通量 MR 筛选矩阵
# ==============================================================================

# 1. 环境初始化
setwd("/Users/bing/DH")

# 检查并安装核心包 (如果还没有安装 TwoSampleMR，系统将自动从 Github 安装)
if (!require("devtools")) install.packages("devtools")
if (!require("TwoSampleMR")) devtools::install_github("MRCIEU/TwoSampleMR")
if (!require("ieugwasr")) install.packages("ieugwasr")

library(TwoSampleMR)
library(ieugwasr)
library(dplyr)

# ⚠️ 请在此处输入您的 IEU OpenGWAS API Token (如果已经在 .Renviron 配置过，可注释掉此行)
# ieugwasr::get_access_token() # 检查token状态

# ==============================================================================
# 2. 定义高通量筛选矩阵 (暴露池 vs 结局池)
# ==============================================================================
cat(">>> 1. 正在初始化 GWAS 数据矩阵池...\n")

# 定义暴露列表 (来源于 NHANES 第一部分发现的睡眠、疲劳、认知等)
exposure_list <- list(
  "Tiredness" = "ukb-b-5262",
  "Insomnia" = "ebi-a-GCST007363",
  "Sleep_duration" = "ebi-a-GCST006460",
  "Major_Depression" = "ebi-a-GCST005902",
  "Cognitive_performance" = "ebi-a-GCST006572"
)

# 定义结局列表 (心血管全家桶)
outcome_list <- list(
  "Heart_Failure" = "ebi-a-GCST009541",
  "Atrial_Fibrillation" = "ebi-a-GCST006414",
  "Coronary_Disease" = "ebi-a-GCST005111",
  "Ischemic_Stroke" = "ebi-a-GCST006908"
)

# 建立一个空的数据框，用于存放全部结果
all_mr_results <- data.frame()

# ==============================================================================
# 3. 嵌套循环执行双样本 MR 筛选 (内建 tryCatch 容错防崩机制)
# ==============================================================================
cat(">>> 2. 开始执行自动化批量 MR 筛选，这可能需要一些时间，请耐心等待...\n")

for (exp_name in names(exposure_list)) {
  exp_id <- exposure_list[[exp_name]]
  cat(sprintf("\n[开始提取暴露数据] %s (ID: %s)...\n", exp_name, exp_id))
  
  # 提取暴露的工具变量 (IVs)，默认 P < 5e-8，并自动进行连锁不平衡剔除 (clumping)
  # 使用 tryCatch 防止某个暴露无显著 SNP 导致循环崩溃
  exp_dat <- tryCatch({
    extract_instruments(outcomes = exp_id)
  }, error = function(e) {
    cat(sprintf("   ⚠️ 暴露 %s 提取失败，跳过。\n", exp_name))
    return(NULL)
  })
  
  if (is.null(exp_dat) || nrow(exp_dat) == 0) {
    cat("   ⚠️ 未找到有效 SNPs，跳过当前暴露。\n")
    next
  }
  
  cat(sprintf("   ✅ 成功提取到 %d 个有效 SNPs。开始匹配结局池...\n", nrow(exp_dat)))
  
  for (out_name in names(outcome_list)) {
    out_id <- outcome_list[[out_name]]
    cat(sprintf("   -> [匹配结局] %s (ID: %s)...\n", out_name, out_id))
    
    tryCatch({
      # 提取结局数据中与 IVs 对应的 SNPs
      out_dat <- extract_outcome_data(
        snps = exp_dat$SNP,
        outcomes = out_id,
        proxies = TRUE # 允许寻找代理 SNP 以防缺失
      )
      
      if (!is.null(out_dat) && nrow(out_dat) > 0) {
        # 协调等位基因 (Harmonise)
        dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat)
        
        # 执行 MR 分析
        res <- mr(dat)
        
        # 提取 IVW (Inverse variance weighted) 方法的核心结果
        ivw_res <- res %>% filter(method == "Inverse variance weighted")
        
        if (nrow(ivw_res) > 0) {
          # 将结果组装并追加到总表
          temp_res <- data.frame(
            Exposure = exp_name,
            Exposure_ID = exp_id,
            Outcome = out_name,
            Outcome_ID = out_id,
            Method = "IVW",
            SNPs_Count = ivw_res$nsnp,
            Beta = ivw_res$b,
            SE = ivw_res$se,
            P_value = ivw_res$pval,
            OR = exp(ivw_res$b),
            OR_LCI95 = exp(ivw_res$b - 1.96 * ivw_res$se),
            OR_UCI95 = exp(ivw_res$b + 1.96 * ivw_res$se)
          )
          all_mr_results <- bind_rows(all_mr_results, temp_res)
          cat(sprintf("      ✔️ 完成! IVW P-value: %g\n", ivw_res$pval))
        }
      } else {
        cat("      ⚠️ 在结局中未找到足够的匹配 SNPs。\n")
      }
    }, error = function(e) {
      cat(sprintf("      ❌ 分析报错，跳过该组合: %s\n", e$message))
    })
  }
}

# ==============================================================================
# 4. 汇总、排序并导出战果
# ==============================================================================
cat("\n>>> 3. 所有筛选完毕！正在导出结果并筛选出潜在的阳性靶点...\n")

if (nrow(all_mr_results) > 0) {
  # 按照 P 值从小到大排序，阳性结果置顶
  all_mr_results <- all_mr_results %>% arrange(P_value)
  
  # 增加一列直观标记阳性
  all_mr_results$Significant <- ifelse(all_mr_results$P_value < 0.05, "Yes (*)", "No")
  
  # 保存为精美的 CSV 筛选总表
  write.csv(all_mr_results, file = "DH_Part2_Step1_MR_Screening_Results.csv", row.names = FALSE)
  
  cat("========================================================\n")
  cat("🎉 第二部分 Step 1 执行完毕！\n")
  cat("您可以在 /Users/bing/DH/ 目录下找到 DH_Part2_Step1_MR_Screening_Results.csv\n")
  cat(sprintf("共完成 %d 组因果关系测试。其中 P < 0.05 的阳性组合有 %d 组！\n", 
              nrow(all_mr_results), sum(all_mr_results$P_value < 0.05)))
  cat("========================================================\n")
} else {
  cat("极其罕见：所有组合均未成功运行，请检查网络或 API Token。\n")
}

# ==============================================================================
# 项目名称: Double Heart (DH)
# 分析模块: 第二部分：双向孟德尔随机化 (MR) 分析
# 步骤名称: Step 2: 核心阳性靶点深度验证与独立单页 PDF 导出 (重构版)
# ==============================================================================

# 1. 环境初始化
setwd("/Users/bing/DH")

# 检查并安装绘图必须的 ggplot2 (如果缺失)
if (!require("ggplot2")) install.packages("ggplot2")

library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(ggplot2)

# ==============================================================================
# 2. 锁定 Step 1 筛选出的“黄金阳性组合”
# ==============================================================================
cat(">>> 1. 正在初始化核心阳性靶点...\n")

# 结局统一为: 心力衰竭 (Heart Failure)
outcome_id <- "ebi-a-GCST009541"
outcome_name <- "Heart_Failure"

# 两个高度显著的暴露
target_exposures <- list(
  "Major_Depression" = "ebi-a-GCST005902",
  "Cognitive_Performance" = "ebi-a-GCST006572"
)

# ==============================================================================
# 3. 自动化高级 MR 分析与独立图形化流水线
# ==============================================================================

for (exp_name in names(target_exposures)) {
  exp_id <- target_exposures[[exp_name]]
  cat(sprintf("\n======================================================\n"))
  cat(sprintf("🔥 正在深度剖析: [%s] ➡️ [%s]\n", exp_name, outcome_name))
  cat(sprintf("======================================================\n"))
  
  # A. 提取与协调 (Harmonise)
  cat("   -> 正在提取工具变量并与结局数据进行协调...\n")
  exp_dat <- extract_instruments(outcomes = exp_id)
  out_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = outcome_id, proxies = TRUE)
  dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat)
  
  # B. 主 MR 分析计算 (提取多种算法结果以确保稳健)
  cat("   -> 正在执行多算法 MR 分析计算...\n")
  res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode"))
  res_or <- generate_odds_ratios(res)
  write.csv(res_or, file = sprintf("DH_Part2_Step2_%s_to_HF_Results.csv", exp_name), row.names = FALSE)
  
  # C. 敏感性分析 (异质性、多效性、留一法)
  cat("   -> 正在进行异质性及水平多效性检验...\n")
  het <- mr_heterogeneity(dat)
  write.csv(het, file = sprintf("DH_Part2_Step2_%s_to_HF_Heterogeneity.csv", exp_name), row.names = FALSE)
  
  pleio <- mr_pleiotropy_test(dat)
  write.csv(pleio, file = sprintf("DH_Part2_Step2_%s_to_HF_Pleiotropy.csv", exp_name), row.names = FALSE)
  
  loo <- mr_leaveoneout(dat)
  
  # D. 渲染并导出独立高清单页 PDF (完美适配 LaTeX 插图排版)
  cat("   -> 正在逐一渲染并导出高阶可视化独立图表...\n")
  
  # 统一主题样式，去背投白底，学术风
  my_theme <- theme_classic() + theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  # 图 1: 散点图 (Scatter)
  pdf(sprintf("DH_Part2_Step2_%s_Scatter.pdf", exp_name), width = 7, height = 5.5)
  p1 <- mr_scatter_plot(res, dat)
  print(p1[[1]] + my_theme + ggtitle(sprintf("Scatter Plot: %s on HF", gsub("_", " ", exp_name))))
  dev.off()
  
  # 图 2: 森林图 (Forest)
  res_single <- mr_singlesnp(dat)
  pdf(sprintf("DH_Part2_Step2_%s_Forest.pdf", exp_name), width = 7, height = 6.5)
  p2 <- mr_forest_plot(res_single)
  print(p2[[1]] + my_theme + ggtitle(sprintf("Forest Plot: %s on HF", gsub("_", " ", exp_name))))
  dev.off()
  
  # 图 3: 留一法 (Leave-one-out)
  pdf(sprintf("DH_Part2_Step2_%s_LOO.pdf", exp_name), width = 7, height = 7)
  p3 <- mr_leaveoneout_plot(loo)
  print(p3[[1]] + my_theme + ggtitle("Leave-one-out Sensitivity Analysis"))
  dev.off()
  
  # 图 4: 漏斗图 (Funnel)
  pdf(sprintf("DH_Part2_Step2_%s_Funnel.pdf", exp_name), width = 7, height = 5.5)
  p4 <- mr_funnel_plot(res_single)
  print(p4[[1]] + my_theme + ggtitle("Funnel Plot for Asymmetry"))
  dev.off()
  
  cat(sprintf("   ✅ 分析完成！四大金刚图已各自独立保存。\n"))
}

cat("\n==============================================================================\n")
cat("SUCCESS! 第二部分 Step 2 图表拆分重建完毕！\n")
cat("所有的图现在都已被单独剥离保存，极其方便您在 LaTeX 中灵活引用。\n")
cat("==============================================================================\n")

# ==============================================================================
# 项目名称: Double Heart (DH)
# 分析模块: 第二部分：双向孟德尔随机化 (MR) 分析
# 步骤名称: Step 3: 反向 MR 分析 (心衰 ➡️ 心理/认知特征) 及独立绘图
# ==============================================================================

# 1. 环境初始化
setwd("/Users/bing/DH")
library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(ggplot2)

# ==============================================================================
# 2. 设定反向 MR 靶点 (暴露与结局互换)
# ==============================================================================
cat(">>> 1. 正在初始化反向 MR 靶点 (Heart Failure 作为暴露)...\n")

# 现在的暴露是: 心力衰竭 (Heart Failure)
exposure_id <- "ebi-a-GCST009541"
exposure_name <- "Heart_Failure"

# 现在的结局是: 重度抑郁与认知能力
target_outcomes <- list(
  "Major_Depression" = "ebi-a-GCST005902",
  "Cognitive_Performance" = "ebi-a-GCST006572"
)

# ==============================================================================
# 3. 提取暴露数据 (Heart Failure)
# ==============================================================================
cat("   -> 正在提取心力衰竭的强效工具变量 (SNPs)...\n")
# 提取心衰 SNPs，默认 P < 5e-8，且自动去除连锁不平衡
exp_dat <- extract_instruments(outcomes = exposure_id)

if (is.null(exp_dat) || nrow(exp_dat) == 0) {
  stop("❌ 提取心衰工具变量失败，请检查网络或 IEU API 状态。")
}
cat(sprintf("   ✅ 成功提取到 %d 个心力衰竭的有效 SNPs。\n", nrow(exp_dat)))

# ==============================================================================
# 4. 循环验证结局，执行反向 MR 分析与可视化
# ==============================================================================

for (out_name in names(target_outcomes)) {
  out_id <- target_outcomes[[out_name]]
  cat(sprintf("\n======================================================\n"))
  cat(sprintf("🔄 正在执行反向验证: [%s] ➡️ [%s]\n", exposure_name, out_name))
  cat(sprintf("======================================================\n"))
  
  # A. 提取结局数据并协调
  cat("   -> 正在提取结局数据并进行 Harmonise 协调...\n")
  out_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = out_id, proxies = TRUE)
  
  if (is.null(out_dat) || nrow(out_dat) == 0) {
    cat(sprintf("   ⚠️ 在 %s 中未找到足够的匹配 SNPs，跳过。\n", out_name))
    next
  }
  
  dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat)
  
  # B. 多算法 MR 分析
  cat("   -> 正在执行多算法 MR 计算...\n")
  res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode"))
  res_or <- generate_odds_ratios(res)
  write.csv(res_or, file = sprintf("DH_Part2_Step3_Reverse_%s_to_%s_Results.csv", exposure_name, out_name), row.names = FALSE)
  
  # C. 敏感性检验
  cat("   -> 正在进行异质性及水平多效性检验...\n")
  het <- mr_heterogeneity(dat)
  write.csv(het, file = sprintf("DH_Part2_Step3_Reverse_%s_to_%s_Heterogeneity.csv", exposure_name, out_name), row.names = FALSE)
  
  pleio <- mr_pleiotropy_test(dat)
  write.csv(pleio, file = sprintf("DH_Part2_Step3_Reverse_%s_to_%s_Pleiotropy.csv", exposure_name, out_name), row.names = FALSE)
  
  loo <- mr_leaveoneout(dat)
  
  # D. 渲染并导出独立高清 PDF (带 Reverse 标记)
  cat("   -> 正在生成反向 MR 独立可视化图表...\n")
  
  my_theme <- theme_classic() + theme(plot.title = element_text(face = "bold", hjust = 0.5))
  clean_out_name <- gsub("_", " ", out_name)
  
  # Scatter Plot
  pdf(sprintf("DH_Part2_Step3_Reverse_%s_Scatter.pdf", out_name), width = 7, height = 5.5)
  p1 <- mr_scatter_plot(res, dat)
  print(p1[[1]] + my_theme + ggtitle(sprintf("Reverse MR: HF on %s", clean_out_name)))
  dev.off()
  
  # Forest Plot
  res_single <- mr_singlesnp(dat)
  pdf(sprintf("DH_Part2_Step3_Reverse_%s_Forest.pdf", out_name), width = 7, height = 6.5)
  p2 <- mr_forest_plot(res_single)
  print(p2[[1]] + my_theme + ggtitle(sprintf("Reverse Forest: HF on %s", clean_out_name)))
  dev.off()
  
  # Leave-one-out
  pdf(sprintf("DH_Part2_Step3_Reverse_%s_LOO.pdf", out_name), width = 7, height = 7)
  p3 <- mr_leaveoneout_plot(loo)
  print(p3[[1]] + my_theme + ggtitle("Reverse Leave-one-out Sensitivity"))
  dev.off()
  
  # Funnel Plot
  pdf(sprintf("DH_Part2_Step3_Reverse_%s_Funnel.pdf", out_name), width = 7, height = 5.5)
  p4 <- mr_funnel_plot(res_single)
  print(p4[[1]] + my_theme + ggtitle("Reverse Funnel Plot"))
  dev.off()
  
  cat(sprintf("   ✅ 反向验证 [%s] 完成！\n", out_name))
}

cat("\n==============================================================================\n")
cat("SUCCESS! 第二部分 Step 3 (反向孟德尔随机化) 圆满完成！\n")
cat("请检查新生成的 Reverse_..._Results.csv 中的 P 值。\n")
cat("==============================================================================\n")

# ==============================================================================
# 项目名称: Double Heart (DH)
# 分析模块: 第二部分：多因素孟德尔随机化 (MVMR) 分析
# 步骤名称: Step 4: 剔除 BMI 混杂因素的独立因果效应验证
# ==============================================================================

# 1. 环境初始化
setwd("/Users/bing/DH")
library(TwoSampleMR)
library(ieugwasr)
library(dplyr)

# ==============================================================================
# 2. 设定 MVMR 的暴露组合与结局
# ==============================================================================
cat(">>> 1. 正在初始化 MVMR 靶点池...\n")

# 结局: 心力衰竭 (Heart Failure)
outcome_id <- "ebi-a-GCST009541"

# 引入经典的混杂因素: BMI (Body mass index, UKB, 460k)
bmi_id <- "ieu-b-40" 

# 需要验证的核心靶点
target_exposures <- list(
  "Major_Depression" = "ebi-a-GCST005902",
  "Cognitive_Performance" = "ebi-a-GCST006572"
)

# 用于保存合并结果的大表
all_mvmr_results <- data.frame()

# ==============================================================================
# 3. 执行多因素 MR 分析 (循环验证两个心理靶点)
# ==============================================================================

for (exp_name in names(target_exposures)) {
  core_exp_id <- target_exposures[[exp_name]]
  
  cat(sprintf("\n======================================================\n"))
  cat(sprintf("🔥 正在执行 MVMR: [%s + BMI] ➡️ [Heart Failure]\n", exp_name))
  cat(sprintf("======================================================\n"))
  
  # A. 提取多个暴露的联合工具变量 (MV Exposures)
  cat("   -> 1/4 正在提取联合暴露 (心理特征 + BMI) 的 SNPs 并去重...\n")
  # mv_extract_exposures 会自动寻找两个暴露的强工具变量，并取并集
  exposure_ids <- c(core_exp_id, bmi_id)
  mv_exp_dat <- mv_extract_exposures(exposure_ids)
  
  if (is.null(mv_exp_dat) || nrow(mv_exp_dat) == 0) {
    cat("   ⚠️ 提取联合暴露 SNPs 失败，跳过。\n")
    next
  }
  
  # B. 提取结局数据
  cat("   -> 2/4 正在结局(Heart Failure)中匹配联合 SNPs...\n")
  mv_out_dat <- extract_outcome_data(snps = mv_exp_dat$SNP, outcomes = outcome_id, proxies = TRUE)
  
  # C. 协调暴露与结局数据 (MV Harmonise)
  cat("   -> 3/4 正在多维度协调等位基因 (MV Harmonise)...\n")
  mv_dat <- mv_harmonise_data(mv_exp_dat, mv_out_dat)
  
  # D. 运行多因素 MR 回归模型
  cat("   -> 4/4 正在计算调整后的独立因果效应...\n")
  res_mvmr <- mv_multiple(mv_dat)
  
  # 提取结果并计算 OR 值
  result_df <- res_mvmr$result %>%
    mutate(
      Model = sprintf("%s_adjusted_for_BMI", exp_name),
      OR = exp(b),
      OR_LCI95 = exp(b - 1.96 * se),
      OR_UCI95 = exp(b + 1.96 * se),
      Significant = ifelse(pval < 0.05, "Yes (*)", "No")
    ) %>%
    select(Model, exposure, outcome, nsnp, b, se, pval, OR, OR_LCI95, OR_UCI95, Significant)
  
  # 将该组的 MVMR 结果追加到总表
  all_mvmr_results <- bind_rows(all_mvmr_results, result_df)
  
  cat(sprintf("   ✅ [%s] 的 MVMR 调整分析完成！\n", exp_name))
}

# ==============================================================================
# 4. 汇总与导出 MVMR 战果
# ==============================================================================
cat("\n>>> 正在整合并导出多因素 MR 最终成果表...\n")

if (nrow(all_mvmr_results) > 0) {
  # 保存为精美的 CSV
  write.csv(all_mvmr_results, file = "DH_Part2_Step4_MVMR_Results.csv", row.names = FALSE)
  
  cat("========================================================\n")
  cat("🎉 伟大的胜利！第二部分 Step 4 (MVMR) 执行完毕！\n")
  cat("您可以在 /Users/bing/DH/ 目录下找到 DH_Part2_Step4_MVMR_Results.csv\n")
  cat("========================================================\n")
} else {
  cat("⚠️ 所有 MVMR 模型均未能生成有效结果，请检查网络或 API。\n")
}
