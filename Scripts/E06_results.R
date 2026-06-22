# ============================================================
# E06_results.R (丰富版)
# 提取扩展队列全部关键数值，包括单变量、多变量、NRI等
# ============================================================

library(dplyr)
library(broom)
library(survival)

ext_dir <- "extended_cohort"
load(file.path(ext_dir, "Analysis_Results_Ext.RData"))
load(file.path(ext_dir, "Clinical_Tools_Ext.RData"))

sink(file.path(ext_dir, "Key_Results_Ext.txt"))

cat("========== 扩展训练集 (NHANES 2009-2018, 5周期) 关键结果 ==========\n\n")
cat("样本量:", nrow(DH_Data), "\n")
cat("心血管死亡事件:", sum(DH_Data$CVD_Death), "\n")
cat("事件率:", round(mean(DH_Data$CVD_Death) * 100, 1), "%\n\n")

# ---- 单变量 Cox ----
cat("--- 单变量 Cox 回归（按 Wald Z 降序） ---\n")
ranking_results <- ranking_results %>% arrange(desc(Z))
for (i in 1:nrow(ranking_results)) {
  cat(sprintf("%-20s  HR=%.3f  Z=%.2f  P=%.4f\n",
              ranking_results$Variable[i],
              ranking_results$HR[i],
              ranking_results$Z[i],
              ranking_results$P[i]))
}
cat("\n")

# ---- 多变量 Cox ----
cat("--- 多变量 Cox 模型 ---\n")
cat("C-index:", round(cox_model_final$concordance["concordance"], 3), "\n")
cat(sprintf("似然比检验: P = %.2e\n\n", summary(cox_model_final)$logtest["pvalue"]))
tidy_cox <- tidy(cox_model_final, exponentiate = TRUE, conf.int = TRUE)
cat("变量                 HR (95% CI)               P\n")
for (i in 1:nrow(tidy_cox)) {
  cat(sprintf("%-15s  %.2f (%.2f-%.2f)           %.4f\n",
              tidy_cox$term[i],
              tidy_cox$estimate[i],
              tidy_cox$conf.low[i],
              tidy_cox$conf.high[i],
              tidy_cox$p.value[i]))
}
cat("\n")

# ---- BMI 非线性检验 ----
if (exists("rcs_fit")) {
  cat("--- BMI 非线性检验 (RCS, 4节点) ---\n")
  anova_rcs <- anova(rcs_fit)
  print(anova_rcs)
  cat(sprintf("BMI 非线性项 P 值: %.4f\n", anova_rcs["BMI", "P"]))
  cat("\n")
}

# ---- PMRS AUC ----
cat("--- PMRS 时间依赖 AUC ---\n")
cat("1-Year AUC:", round(roc_1y$AUC[2], 3), "\n")
cat("3-Year AUC:", round(roc_3y$AUC[2], 3), "\n")
cat("5-Year AUC:", round(roc_5y$AUC[2], 3), "\n")
cat("连续风险评分 3-Year AUC:", round(roc_risk_3y$AUC[2], 3), "\n\n")

# ---- NRI ----
cat("--- 净重分类改善 (NRI) at 36 months ---\n")
if (exists("nri") && !is.null(nri) && "nri" %in% names(nri)) {
  nri_mat <- nri$nri
  extract_val <- function(row_name, col_name) {
    if (row_name %in% rownames(nri_mat) && col_name %in% colnames(nri_mat)) {
      return(round(nri_mat[row_name, col_name], 3))
    } else { return(NA) }
  }
  cat("NRI (总体):", extract_val("NRI", "Estimate"),
      " [95% CI:", extract_val("NRI", "Lower"), ",", extract_val("NRI", "Upper"), "]\n")
  cat("NRI+ (事件):", extract_val("NRI+", "Estimate"), "\n")
  cat("NRI- (无事件):", extract_val("NRI-", "Estimate"), "\n")
  cat("Pr(Up|Case):", extract_val("Pr(Up|Case)", "Estimate"), "\n")
  cat("Pr(Down|Ctrl):", extract_val("Pr(Down|Ctrl)", "Estimate"), "\n")
} else {
  cat("NRI 结果未计算。\n")
}
cat("\n")

# ---- 评分表变量 ----
cat("--- PMRS 评分表变量 ---\n")
cat("Age, Gender, BMI, Cognitive_Score, DPQ030, DPQ040\n")
cat("完整评分表见 output/Table_PMRS_Score.csv\n")

sink()
cat("结果已保存至 extended_cohort/Key_Results_Ext.txt\n")
