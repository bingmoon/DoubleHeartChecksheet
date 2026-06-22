# ============================================================
# E04_clinical_tools.R
# Nomogram, PMRS, DCA, NRI, 时间依赖ROC
# 输入：extended_cohort/Analysis_Results_Ext.RData
# 输出：extended_cohort/Clinical_Tools_Ext.RData
# ============================================================

library(rms)
library(timeROC)
library(nricens)
library(dplyr)
library(survival)

ext_dir <- "extended_cohort"
load(file.path(ext_dir, "Analysis_Results_Ext.RData"))

# PMRS 计算
calc_pmrs <- function(data) {
  score <- integer(nrow(data))
  score <- score + findInterval(data$Age, c(50, 60, 70))
  score <- score + ifelse(data$Gender == "Male", 2, 0)
  score <- score - findInterval(data$BMI, c(25, 30, 35))
  score <- score + findInterval(data$Cognitive_Score, c(4, 8))
  score <- score + data$DPQ030
  dpq040_pts <- c(0, 1, 1, 2)
  score <- score + dpq040_pts[data$DPQ040 + 1]
  return(score)
}

DH_Data$PMRS <- calc_pmrs(DH_Data)

# Nomogram # 核心修复：在这里强行给变量贴上要在列线图上显示的标签 
label(DH_Data$Cognitive_Score) <- "DCS" 
label(DH_Data$DPQ030) <- "Sleep Disturbance" 
label(DH_Data$DPQ040) <- "Fatigue" 

dd <- datadist(DH_Data); options(datadist = "dd") 
cph_model <- cph(cox_formula_final, data = DH_Data, surv = TRUE, x = TRUE, y = TRUE) 
nom <- nomogram(cph_model, fun = list(function(x) x), funlabel = "3-Year Survival", fun.at = c(0.5, 0.7, 0.9), lp = FALSE, conf.int = FALSE)

# 时间依赖ROC
roc_1y <- timeROC(T = DH_Data$Survival_Months, delta = DH_Data$CVD_Death,
                  marker = DH_Data$PMRS, cause = 1, times = 12)
roc_3y <- timeROC(T = DH_Data$Survival_Months, delta = DH_Data$CVD_Death,
                  marker = DH_Data$PMRS, cause = 1, times = 36)
roc_5y <- timeROC(T = DH_Data$Survival_Months, delta = DH_Data$CVD_Death,
                  marker = DH_Data$PMRS, cause = 1, times = 60)
roc_risk_3y <- timeROC(T = DH_Data$Survival_Months, delta = DH_Data$CVD_Death,
                       marker = DH_Data$risk_score, cause = 1, times = 36)

cat("PMRS AUC: 1yr =", round(roc_1y$AUC[2], 3),
    "3yr =", round(roc_3y$AUC[2], 3),
    "5yr =", round(roc_5y$AUC[2], 3), "\n")

# DCA
DH_Data$event_36m <- ifelse(DH_Data$Survival_Months <= 36 & DH_Data$CVD_Death == 1, 1, 0)
risk_death_full <- 1 - predict(cox_model_final, type = "survival", times = 36)
cox_base <- coxph(Surv(Survival_Months, CVD_Death) ~ Age + Gender + BMI, data = DH_Data)
risk_death_base <- 1 - predict(cox_base, type = "survival", times = 36)

thresholds <- seq(0.01, 0.50, by = 0.01)
calc_nb <- function(th, pred, event) {
  pred_bin <- ifelse(pred >= th, 1, 0)
  tp <- sum(pred_bin == 1 & event == 1)
  fp <- sum(pred_bin == 1 & event == 0)
  n <- length(event)
  (tp/n) - (fp/n) * (th / (1 - th))
}

nb_treat_all <- sapply(thresholds, function(th) mean(DH_Data$event_36m) - (1-mean(DH_Data$event_36m)) * th/(1-th))
nb_full <- sapply(thresholds, function(th) calc_nb(th, risk_death_full, DH_Data$event_36m))
nb_base <- sapply(thresholds, function(th) calc_nb(th, risk_death_base, DH_Data$event_36m))

df_dca <- rbind(
  data.frame(threshold = thresholds, net_benefit = nb_treat_all, Model = "Treat All"),
  data.frame(threshold = thresholds, net_benefit = 0, Model = "Treat None"),
  data.frame(threshold = thresholds, net_benefit = nb_full, Model = "Full Model"),
  data.frame(threshold = thresholds, net_benefit = nb_base, Model = "Base Model")
)

# NRI
nri <- nricens(time = DH_Data$Survival_Months, event = DH_Data$CVD_Death,
               p.std = risk_death_base, p.new = risk_death_full,
               t0 = 36, cut = c(0.05, 0.10), niter = 1000)

save(DH_Data, nom, roc_1y, roc_3y, roc_5y, roc_risk_3y, df_dca, nri,
     file = file.path(ext_dir, "Clinical_Tools_Ext.RData"))
message("E04_clinical_tools.R completed.")

# 补丁：生成最新 PMRS 评分表与 NRI 文件
load("extended_cohort/Clinical_Tools_Ext.RData")

# ---- 重新生成 PMRS 评分表 (与 E04 中 make_score_table 完全一致) ----
make_score_table <- function() {
  tab <- rbind(
    data.frame(Variable = "Age",             Category = c("40-49","50-59","60-69","70+"), Points = c(0,1,2,3)),
    data.frame(Variable = "Gender",          Category = c("Female","Male"),               Points = c(0,2)),
    data.frame(Variable = "BMI",             Category = c("<25","25-30","30-35",">35"),    Points = c(0,-1,-2,-3)),
    data.frame(Variable = "Cognitive_Score", Category = c("0-3","4-7","8+"),               Points = c(0,1,2)),
    data.frame(Variable = "DPQ030",          Category = c("0","1","2","3"),                Points = c(0,1,2,3)),
    data.frame(Variable = "DPQ040",          Category = c("0","1","2","3"),                Points = c(0,1,1,2))
  )
  return(tab)
}
PMRS_Score_Table <- make_score_table()

# 保存评分表 CSV
write.csv(PMRS_Score_Table, "extended_cohort/Table_PMRS_Score.csv", row.names = FALSE)
cat("✅ 已生成 Table_PMRS_Score.csv\n")

# ---- 输出 NRI 摘要 ----
sink("extended_cohort/Table_NRI_Summary.txt")
cat("NRI analysis at 36 months (Extended Cohort, N = 2,308, 363 events)\n\n")
if (exists("nri") && !is.null(nri) && "nri" %in% names(nri)) {
  print(nri)
} else {
  cat("NRI object not found.\n")
}
sink()
cat("✅ 已生成 Table_NRI_Summary.txt\n")
