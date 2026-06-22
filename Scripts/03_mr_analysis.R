# ============================================================
# 03_mr_analysis.R
# 双向孟德尔随机化 + 多变量 MR (心理→心衰，代谢中介)
# 输入：环境变量 OPENGWAS_JWT 已设置
# 输出：output/MR_Results.RData, output/MR_Summary.csv, 图表
# ============================================================

library(TwoSampleMR)
library(dplyr)
library(ggplot2)

# ---- 0. Token 检查 ----
if (Sys.getenv("OPENGWAS_JWT") == "") {
  stop("请先设置 OPENGWAS_JWT 环境变量！运行: Sys.setenv(OPENGWAS_JWT = '你的token')")
}
cat("Token 环境变量已就绪。\n")

# ---- 1. 定义 GWAS ID ----
cat("\n>>> 1. 定义暴露与结局矩阵...\n")

exposure_list <- list(
  "Tiredness"             = "ukb-b-5262",
  "Insomnia"              = "ebi-a-GCST007363",
  "Sleep_duration"        = "ebi-a-GCST006460",
  "Major_Depression"      = "ebi-a-GCST005902",
  "Cognitive_performance" = "ebi-a-GCST006572"
)

outcome_list <- list(
  "Heart_Failure"       = "ebi-a-GCST009541",
  "Atrial_Fibrillation" = "ebi-a-GCST006414",
  "Coronary_Disease"    = "ebi-a-GCST005111",
  "Ischemic_Stroke"     = "ebi-a-GCST006908"
)

hf_outcome_id <- "ebi-a-GCST009541"
bmi_id <- "ieu-b-40"

# ---- 2. 高通量正向 MR 筛选矩阵 ----
cat("\n>>> 2. 执行暴露-结局全矩阵 MR 筛选...\n")

all_mr_results <- data.frame()

for (exp_name in names(exposure_list)) {
  exp_id <- exposure_list[[exp_name]]
  cat(sprintf("\n--- 暴露：%s (ID: %s) ---\n", exp_name, exp_id))
  
  exp_dat <- tryCatch({
    extract_instruments(outcomes = exp_id)
  }, error = function(e) {
    cat("  提取暴露 IV 失败：", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(exp_dat) || nrow(exp_dat) == 0) {
    cat("  未提取到有效工具变量，跳过。\n")
    next
  }
  cat(sprintf("  提取到 %d 个工具变量\n", nrow(exp_dat)))
  
  for (out_name in names(outcome_list)) {
    out_id <- outcome_list[[out_name]]
    cat(sprintf("  -> 结局：%s (ID: %s)\n", out_name, out_id))
    
    tryCatch({
      out_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = out_id, proxies = TRUE)
      if (is.null(out_dat) || nrow(out_dat) == 0) next
      
      dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat)
      res <- mr(dat)
      
      ivw_res <- res %>% filter(method == "Inverse variance weighted")
      if (nrow(ivw_res) > 0) {
        all_mr_results <<- bind_rows(all_mr_results, data.frame(
          Exposure  = exp_name,
          Outcome   = out_name,
          SNP_count = ivw_res$nsnp,
          Beta      = ivw_res$b,
          SE        = ivw_res$se,
          P_value   = ivw_res$pval,
          OR        = exp(ivw_res$b),
          OR_lower  = exp(ivw_res$b - 1.96 * ivw_res$se),
          OR_upper  = exp(ivw_res$b + 1.96 * ivw_res$se)
        ))
        cat(sprintf("    IVW p = %.4f, OR = %.2f\n", ivw_res$pval, exp(ivw_res$b)))
      }
    }, error = function(e) {
      cat("    分析报错：", e$message, "\n")
    })
  }
}

if (nrow(all_mr_results) > 0) {
  all_mr_results <- all_mr_results %>% arrange(P_value)
  write.csv(all_mr_results, "output/MR_Screening_Matrix.csv", row.names = FALSE)
  print(head(all_mr_results, 10))
}

# ---- 3. 核心靶点深度验证：Major Depression / Cognitive Performance → Heart Failure ----
cat("\n>>> 3. 核心阳性靶点深度验证...\n")

target_exposures <- list(
  "Major_Depression"      = "ebi-a-GCST005902",
  "Cognitive_Performance" = "ebi-a-GCST006572"
)

mr_forward_list <- list()

for (exp_name in names(target_exposures)) {
  exp_id <- target_exposures[[exp_name]]
  cat(sprintf("\n===== 正向 MR：%s → Heart Failure =====\n", exp_name))
  
  exp_dat <- extract_instruments(outcomes = exp_id)
  out_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = hf_outcome_id, proxies = TRUE)
  dat <- harmonise_data(exp_dat, out_dat)
  
  res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", 
                                 "mr_weighted_median", "mr_simple_mode"))
  res_or <- generate_odds_ratios(res)
  het   <- mr_heterogeneity(dat)
  pleio <- mr_pleiotropy_test(dat)
  loo   <- mr_leaveoneout(dat)
  res_single <- mr_singlesnp(dat)
  
  mr_forward_list[[exp_name]] <- list(
    dat = dat, res = res, res_or = res_or,
    het = het, pleio = pleio, loo = loo, res_single = res_single
  )
  
  write.csv(res_or, paste0("output/MR_", exp_name, "_to_HF.csv"), row.names = FALSE)
  write.csv(het,   paste0("output/MR_", exp_name, "_heterogeneity.csv"), row.names = FALSE)
  write.csv(pleio, paste0("output/MR_", exp_name, "_pleiotropy.csv"), row.names = FALSE)
  
  my_theme <- theme_classic() + theme(plot.title = element_text(face = "bold", hjust = 0.5))
  clean_name <- gsub("_", " ", exp_name)
  
  pdf(paste0("output/MR_", exp_name, "_Scatter.pdf"), width = 7, height = 5.5)
  p1 <- mr_scatter_plot(res, dat)
  print(p1[[1]] + my_theme + ggtitle(paste0("Scatter: ", clean_name, " on Heart Failure")))
  dev.off()
  
  pdf(paste0("output/MR_", exp_name, "_Forest.pdf"), width = 7, height = 6.5)
  p2 <- mr_forest_plot(res_single)
  print(p2[[1]] + my_theme + ggtitle(paste0("Forest: ", clean_name, " SNPs")))
  dev.off()
  
  pdf(paste0("output/MR_", exp_name, "_LOO.pdf"), width = 7, height = 7)
  p3 <- mr_leaveoneout_plot(loo)
  print(p3[[1]] + my_theme + ggtitle("Leave-one-out Sensitivity"))
  dev.off()
  
  pdf(paste0("output/MR_", exp_name, "_Funnel.pdf"), width = 7, height = 5.5)
  p4 <- mr_funnel_plot(res_single)
  print(p4[[1]] + my_theme + ggtitle("Funnel Plot"))
  dev.off()
  
  cat(sprintf("  %s 验证完成，图表已保存。\n", exp_name))
}

# ---- 4. 反向 MR：Heart Failure → 心理特质 ----
cat("\n>>> 4. 反向 MR 验证...\n")

exp_dat_hf <- extract_instruments(outcomes = hf_outcome_id)
mr_reverse_list <- list()

if (!is.null(exp_dat_hf) && nrow(exp_dat_hf) > 0) {
  for (out_name in names(target_exposures)) {
    out_id <- target_exposures[[out_name]]
    cat(sprintf("\n反向：Heart Failure → %s\n", out_name))
    
    tryCatch({
      out_dat <- extract_outcome_data(snps = exp_dat_hf$SNP, outcomes = out_id, proxies = TRUE)
      if (is.null(out_dat) || nrow(out_dat) == 0) next
      
      dat <- harmonise_data(exp_dat_hf, out_dat)
      res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                                     "mr_weighted_median", "mr_simple_mode"))
      res_or <- generate_odds_ratios(res)
      
      mr_reverse_list[[out_name]] <- list(dat = dat, res = res, res_or = res_or)
      write.csv(res_or, paste0("output/MR_Reverse_HF_to_", out_name, ".csv"), row.names = FALSE)
      
      pdf(paste0("output/MR_Reverse_HF_to_", out_name, "_Scatter.pdf"), width = 7, height = 5.5)
      p <- mr_scatter_plot(res, dat)
      print(p[[1]] + my_theme + ggtitle(paste("Reverse: HF on", gsub("_", " ", out_name))))
      dev.off()
    }, error = function(e) {
      cat("  反向 MR 失败：", e$message, "\n")
    })
  }
}

# ---- 5. 多变量 MR (MVMR)：Major Depression + BMI → Heart Failure ----
cat("\n>>> 5. MVMR 分析...\n")

tryCatch({
  mv_exposure_ids <- c(target_exposures[["Major_Depression"]], bmi_id)
  mv_exp <- mv_extract_exposures(mv_exposure_ids)
  mv_out <- extract_outcome_data(snps = mv_exp$SNP, outcomes = hf_outcome_id, proxies = TRUE)
  mv_dat <- mv_harmonise_data(mv_exp, mv_out)
  res_mvmr <- mv_multiple(mv_dat)
  
  mvmr_table <- res_mvmr$result %>%
    mutate(OR = exp(b), OR_lower = exp(b - 1.96*se), OR_upper = exp(b + 1.96*se))
  write.csv(mvmr_table, "output/MVMR_Depression_BMI_HF.csv", row.names = FALSE)
  print(mvmr_table)
}, error = function(e) {
  cat("MVMR 失败：", e$message, "\n")
  res_mvmr <- NULL
})

# ---- 6. 保存全部结果 ----
save(all_mr_results, mr_forward_list, mr_reverse_list, res_mvmr,
     file = "output/MR_Results.RData")

cat("\n============================================================\n")
cat("SUCCESS！03_mr_analysis.R 执行完毕。\n")
cat("============================================================\n")

# ---- 7. 汇总 MR 结果 ----
load("output/MR_Results.RData")

# 筛选矩阵中针对 Heart Failure 的结果
hf_results <- all_mr_results %>%
  filter(Outcome == "Heart_Failure") %>%
  mutate(Source = "MR Screening (IVW)")

# 核心靶点深度验证结果（取 IVW）
forward_summary <- lapply(names(mr_forward_list), function(exp) {
  res <- mr_forward_list[[exp]]$res
  ivw <- res %>% filter(method == "Inverse variance weighted")
  data.frame(
    Exposure = exp,
    Outcome = "Heart Failure",
    Source = "Forward MR (IVW)",
    SNP_count = ivw$nsnp,
    Beta = ivw$b,
    SE = ivw$se,
    P_value = ivw$pval,
    OR = exp(ivw$b),
    OR_lower = exp(ivw$b - 1.96 * ivw$se),
    OR_upper = exp(ivw$b + 1.96 * ivw$se)
  )
}) %>% bind_rows()

# 反向 MR 结果（取 IVW）
reverse_summary <- lapply(names(mr_reverse_list), function(out) {
  res <- mr_reverse_list[[out]]$res
  ivw <- res %>% filter(method == "Inverse variance weighted")
  data.frame(
    Exposure = "Heart Failure",
    Outcome = out,
    Source = "Reverse MR (IVW)",
    SNP_count = ivw$nsnp,
    Beta = ivw$b,
    SE = ivw$se,
    P_value = ivw$pval,
    OR = exp(ivw$b),
    OR_lower = exp(ivw$b - 1.96 * ivw$se),
    OR_upper = exp(ivw$b + 1.96 * ivw$se)
  )
}) %>% bind_rows()

# MVMR 结果
if (!is.null(res_mvmr)) {
  mvmr_summary <- res_mvmr$result %>%
    mutate(
      Exposure = exposure,
      Outcome = "Heart Failure",
      Source = "MVMR (adjusted for BMI)",
      OR = exp(b),
      OR_lower = exp(b - 1.96 * se),
      OR_upper = exp(b + 1.96 * se)
    ) %>%
    select(Exposure, Outcome, Source, nsnp, b, se, pval, OR, OR_lower, OR_upper) %>%
    rename(Beta = b, SE = se, P_value = pval, SNP_count = nsnp)
} else {
  mvmr_summary <- NULL
}

# 合并并保存
mr_summary <- bind_rows(hf_results, forward_summary, reverse_summary, mvmr_summary)

write.csv(mr_summary, "output/MR_Summary.csv", row.names = FALSE)
print(mr_summary, digits = 4)
cat("\n汇总完成，已保存至 output/MR_Summary.csv\n")
