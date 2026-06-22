# ============================================================
# E02_analysis.R
# 基于扩展队列的单变量 + 多变量 Cox 模型 + BMI 非线性检验
# 输入：extended_cohort/Extended_Cohort_Data.RData
# 输出：extended_cohort/Analysis_Results_Ext.RData
# ============================================================

library(survival)
library(dplyr)
library(rms)

ext_dir <- "extended_cohort"
load(file.path(ext_dir, "Extended_Cohort_Data.RData"))

DH_Data <- DH_Data_Ext
cat("Cohort size:", nrow(DH_Data), "Events:", sum(DH_Data$CVD_Death), "\n")

# ---- 单变量筛选 ----
psych_vars <- c("PHQ9_Total", "Cognitive_Score", "DPQ070", "DPQ030", "DPQ040")
ranking_list <- list()

for (var in psych_vars) {
  if (all(is.na(DH_Data[[var]]))) next
  formula <- as.formula(paste("Surv(Survival_Months, CVD_Death) ~", var))
  fit <- try(coxph(formula, data = DH_Data), silent = TRUE)
  if (inherits(fit, "try-error")) next
  s <- summary(fit)$coefficients
  ranking_list[[var]] <- data.frame(
    Variable = var, HR = exp(s[1,"coef"]), Z = abs(s[1,"z"]), P = s[1,"Pr(>|z|)"]
  )
}

ranking_results <- bind_rows(ranking_list) %>% arrange(desc(Z))
print(ranking_results)

# ---- 多变量 Cox 模型 ----
base_vars <- "Age + Gender + BMI"
psych_vars_model <- c("Cognitive_Score", "DPQ030", "DPQ040")

formula_str <- paste("Surv(Survival_Months, CVD_Death) ~", base_vars, "+",
                     paste(psych_vars_model, collapse = " + "))
cox_formula_final <- as.formula(formula_str)
cox_model_final <- coxph(cox_formula_final, data = DH_Data)
print(summary(cox_model_final))
cat("C-index:", cox_model_final$concordance["concordance"], "\n")

# ---- 风险评分与分层 ----
DH_Data$risk_score <- predict(cox_model_final, type = "lp")
DH_Data$risk_group <- ifelse(DH_Data$risk_score >= median(DH_Data$risk_score), "High", "Low")
km_fit_final <- survfit(Surv(Survival_Months, CVD_Death) ~ risk_group, data = DH_Data)

# ---- BMI 非线性检验（限制性立方样条） ----
cat("\n========== BMI 非线性检验 (RCS) ==========\n") 
dd <- datadist(DH_Data) 
options(datadist = "dd")

# 使用 rms::cph 拟合含 RCS 的模型，4个节点
rcs_fit <- cph(Surv(Survival_Months, CVD_Death) ~ rcs(BMI, 4) + Age + Gender + DPQ040,
               data = DH_Data, surv = TRUE)
# 输出方差分析（含非线性项检验）
print(anova(rcs_fit))
# 提取非线性项 P 值
nonlin_p <- anova(rcs_fit)["BMI", "P"]
cat(sprintf("BMI 非线性项 P 值: %.4f\n", nonlin_p))
if (nonlin_p > 0.05) {
  cat("结论：BMI 与心血管死亡的关联无显著非线性趋势（P > 0.05），不支持 U 型曲线假设。\n")
} else {
  cat("结论：BMI 与心血管死亡存在显著非线性关联，需进一步描述曲线形态。\n")
}

# ---- 保存 ----
save(ranking_results, cox_model_final, cox_formula_final, km_fit_final, DH_Data,
     psych_vars, psych_vars_model, rcs_fit, file = file.path(ext_dir, "Analysis_Results_Ext.RData"))
message("E02_analysis.R completed.")
