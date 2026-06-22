# ============================================================
# E05_figures.R
# 生成所有图表（扩展队列观察性 + MR + BMI非线性曲线）
# 前提：已运行 E01, E02, E03 (MR), E04
# 输出：extended_cohort/Figures/
# ============================================================

library(ggplot2)
library(survminer)
library(dplyr)
library(broom)
library(rms)
library(survival)
library(TwoSampleMR)

ext_dir <- "extended_cohort"
fig_dir <- file.path(ext_dir, "Figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ---- 统一 Nature 风格主题 ----
theme_nature <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

# ---- 加载观察性分析结果 ----
load(file.path(ext_dir, "Analysis_Results_Ext.RData"))
load(file.path(ext_dir, "Clinical_Tools_Ext.RData"))

# ==============================================================================
# 第一部分：观察性分析图表
# ==============================================================================

# 统一变量名 DCS
ranking_results$Variable <- gsub("Cognitive_Score", "DCS", ranking_results$Variable)

# ---- Fig1: Lollipop (样本量已更新为2308) ----
p1 <- ggplot(ranking_results, aes(x = reorder(Variable, Z), y = Z)) +
  geom_segment(aes(xend = Variable, y = 0), color = "grey60", linewidth = 1) +
  geom_point(size = 5, color = "#0072B5") +
  geom_text(aes(label = sprintf("Z=%.2f", Z)), hjust = -0.3, size = 3.5) +
  coord_flip(ylim = c(0, max(ranking_results$Z) * 1.3)) +
  labs(title = "Prognostic Signals (Extended Cohort)",
       subtitle = "NHANES 2009-2018, N=2,308, 363 events",
       x = NULL, y = "Absolute Wald Z-score") + theme_nature
ggsave(file.path(fig_dir, "Fig1_Lollipop.pdf"), p1, width = 8, height = 5)

# ---- Fig2: Forest ----
tidy_cox <- tidy(cox_model_final, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)")
tidy_cox$term <- gsub("Cognitive_Score", "DCS", tidy_cox$term)

p2 <- ggplot(tidy_cox, aes(x = estimate, y = reorder(term, estimate))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 4, color = "#0072B5") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.25, linewidth = 1, color = "#0072B5") +
  scale_x_log10(breaks = c(0.5, 1, 2)) +
  labs(title = "Multivariable Predictors (Extended Cohort)",
       subtitle = paste("C-index =", round(cox_model_final$concordance["concordance"], 3)),
       x = "Hazard Ratio (95% CI)", y = NULL) + theme_nature
ggsave(file.path(fig_dir, "Fig2_Forest.pdf"), p2, width = 8, height = 6)

# ---- Fig3: KM ----
km_plot <- ggsurvplot(km_fit_final, data = DH_Data, pval = TRUE, 
palette = c("#0072B5", "#BC3C29"), 
legend.title = "Risk Group", 
legend.labs = c("Low Risk", "High Risk"), 
xlab = "Survival Months", ylab = "Survival Probability", 
censor = FALSE, 
title = "Kaplan-Meier (Extended Cohort)")

km_plot$plot <- km_plot$plot + theme_nature
ggsave(file.path(fig_dir, "Fig3_KM.pdf"), km_plot$plot, width = 8, height = 6)

# ---- Fig4: DCA ----
p4 <- ggplot(df_dca, aes(x = threshold, y = net_benefit, color = Model, linetype = Model)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("Treat All" = "#E18727", "Treat None" = "#20854E",
                                "Full Model" = "#0072B5", "Base Model" = "#BC3C29")) +
  scale_linetype_manual(values = c("Treat All" = "dashed", "Treat None" = "dashed",
                                   "Full Model" = "solid", "Base Model" = "solid")) +
  coord_cartesian(xlim = c(0, 0.3), ylim = c(-0.02, 0.15)) +
  labs(title = "DCA at 36 Months (Extended Cohort)",
       x = "Threshold Probability", y = "Net Benefit") + theme_nature
ggsave(file.path(fig_dir, "Fig4_DCA.pdf"), p4, width = 8, height = 6)

# ---- Fig5: TimeROC ----
pdf(file.path(fig_dir, "Fig5_TimeROC.pdf"), width = 7, height = 7)
plot(roc_1y, time = 12, col = "#E18727", lwd = 2.5, add = FALSE,
     main = "Time-Dependent ROC for PMRS (Extended Cohort)",
     xlab = "1 - Specificity", ylab = "Sensitivity")
plot(roc_3y, time = 36, col = "#0072B5", lwd = 2.5, add = TRUE)
plot(roc_5y, time = 60, col = "#BC3C29", lwd = 2.5, add = TRUE)
abline(0, 1, lty = 2)
legend("bottomright",
       legend = c(paste0("1-Year AUC=", round(roc_1y$AUC[2], 3)),
                  paste0("3-Year AUC=", round(roc_3y$AUC[2], 3)),
                  paste0("5-Year AUC=", round(roc_5y$AUC[2], 3))),
       col = c("#E18727", "#0072B5", "#BC3C29"), lwd = 2.5, bty = "n")
dev.off()

# ---- Fig6: Nomogram ----
pdf(file.path(fig_dir, "Fig6_Nomogram.pdf"), width = 14, height = 9)
plot(nom, cex.axis = 0.6, cex.var = 0.7, lmgp = 0.2, xfrac = 0.35, col.grid = gray(0.9))
title(main = "Nomogram for 3-Year Cardiovascular Survival (Extended Cohort)", cex.main = 1.2)
dev.off()

# ---- 校准曲线 ----
dd <- datadist(DH_Data); options(datadist = "dd")
fit_rms <- cph(cox_formula_final, data = DH_Data, surv = TRUE, x = TRUE, y = TRUE)
pdf(file.path(fig_dir, "Fig5_Calibration.pdf"), width = 7, height = 7)
cal <- calibrate(fit_rms, method = "boot", B = 1000, u = 36)
plot(cal, xlab = "Predicted 3-Year Survival", ylab = "Observed 3-Year Survival",
     main = "Calibration Curve (Extended Cohort)")
dev.off()

# ---- PMRS KM ----
DH_Data$PMRS_Risk <- ifelse(DH_Data$PMRS >= 8, "High Risk (>=8)", "Low Risk (<8)")
km_pmrs <- survfit(Surv(Survival_Months, CVD_Death) ~ PMRS_Risk, data = DH_Data)

km_pmrs_plot <- ggsurvplot(km_pmrs, data = DH_Data, pval = TRUE, 
palette = c("#0072B5", "#BC3C29"), 
legend.title = "PMRS", 
xlab = "Survival Months", ylab = "Survival Probability", 
censor = FALSE, # 核心修复：关闭删失点显示 
title = "Kaplan-Meier by PMRS (Extended Cohort)")

ggsave(file.path(fig_dir, "Fig3_KM_Checksheet.pdf"), km_pmrs_plot$plot, width = 8, height = 6)

# ---- BMI 非线性曲线 (RCS) ----
# 提取预测值
bmi_range <- quantile(DH_Data$BMI, probs = c(0.01, 0.99), na.rm = TRUE)
bmi_seq <- seq(bmi_range[1], bmi_range[2], length.out = 100)
pred_data <- data.frame(BMI = bmi_seq, Age = median(DH_Data$Age), Gender = "Male", DPQ040 = 0)
pred <- Predict(rcs_fit, BMI = bmi_seq, Age = median(DH_Data$Age), Gender = "Male", DPQ040 = 0, fun = exp, conf.int = TRUE)
pred_df <- as.data.frame(pred)

p_rcs <- ggplot(pred_df, aes(x = BMI, y = yhat)) + geom_line(color = "#0072B5", linewidth = 1.2) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#0072B5", alpha = 0.15) + geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") + labs(title = "Restricted Cubic Spline: BMI and CV Mortality",  
subtitle = sprintf("Nonlinear P = %.4f", anova(rcs_fit)[" Nonlinear", "P"]), x = "Body Mass Index (kg/m²)", y = "Hazard Ratio (vs. BMI = 25)") + theme_nature
ggsave(file.path(fig_dir, "FigS7_BMI_RCS.pdf"), p_rcs, width = 8, height = 5)

# ==============================================================================
# 第二部分：孟德尔随机化图表（完全自包含生成）
# ==============================================================================

if (!file.exists("output/MR_Results.RData")) {
  stop("MR结果文件 output/MR_Results.RData 不存在。请先运行 03_mr_analysis.R。")
}
load("output/MR_Results.RData")

plot_mr_scatter <- function(dat, res, title) {
  p <- mr_scatter_plot(res, dat)
  p[[1]] + theme_nature + ggtitle(title) + theme(legend.position = "none")
}

plot_mr_forest <- function(res_single, title) {
  p <- mr_forest_plot(res_single)
  p[[1]] + theme_nature + ggtitle(title) + theme(legend.position = "none")
}

plot_mr_loo <- function(loo, title) {
  p <- mr_leaveoneout_plot(loo)
  p[[1]] + theme_nature + ggtitle(title) + theme(legend.position = "none")
}

plot_mr_funnel <- function(res_single, title) {
  p <- mr_funnel_plot(res_single)
  p[[1]] + theme_nature + ggtitle(title) + theme(legend.position = "none")
}

# ---- Major Depression → Heart Failure ----
if (!is.null(mr_forward_list[["Major_Depression"]])) {
  md <- mr_forward_list[["Major_Depression"]]
  ggsave(file.path(fig_dir, "Fig7_MDD_Scatter.pdf"),
         plot_mr_scatter(md$dat, md$res, "Major Depression on Heart Failure"), width = 7, height = 5.5)
  ggsave(file.path(fig_dir, "Fig8_MDD_Forest.pdf"),
         plot_mr_forest(md$res_single, "Major Depression SNPs on HF"), width = 7, height = 6.5)
  ggsave(file.path(fig_dir, "FigS1_MDD_LOO.pdf"),
         plot_mr_loo(md$loo, "Leave-one-out: MDD on HF"), width = 7, height = 7)
  ggsave(file.path(fig_dir, "FigS2_MDD_Funnel.pdf"),
         plot_mr_funnel(md$res_single, "Funnel Plot: MDD on HF"), width = 7, height = 5.5)
}

# ---- Cognitive Performance → Heart Failure ----
if (!is.null(mr_forward_list[["Cognitive_Performance"]])) { cp <- mr_forward_list[["Cognitive_Performance"]] 

ggsave(file.path(fig_dir, "Fig9_CogPerf_Scatter.pdf"), plot_mr_scatter(cp$dat, cp$res, "Cognitive Performance on Heart Failure"), width = 7, height = 5.5) # 
 
ggsave(file.path(fig_dir, "Fig10_CogPerf_Forest.pdf"), plot_mr_forest(cp$res_single, "Cognitive Performance SNPs on HF"), width = 7, height = 20) 
 
ggsave(file.path(fig_dir, "FigS3_CogPerf_LOO.pdf"), plot_mr_loo(cp$loo, "Leave-one-out: CogPerf on HF"), width = 7, height = 20) 

ggsave(file.path(fig_dir, "FigS4_CogPerf_Funnel.pdf"), plot_mr_funnel(cp$res_single, "Funnel Plot: CogPerf on HF"), width = 7, height = 5.5) 
}

# ---- 反向 MR ----
if (!is.null(mr_reverse_list[["Major_Depression"]])) {
  rev_md <- mr_reverse_list[["Major_Depression"]]
  p_rev <- mr_scatter_plot(rev_md$res, rev_md$dat)
  ggsave(file.path(fig_dir, "FigS5_Reverse_MDD_Scatter.pdf"),
         p_rev[[1]] + theme_nature + ggtitle("Reverse MR: HF on Major Depression"), width = 7, height = 6.5)
}

if (!is.null(mr_reverse_list[["Cognitive_Performance"]])) {
  rev_cp <- mr_reverse_list[["Cognitive_Performance"]]
  p_rev2 <- mr_scatter_plot(rev_cp$res, rev_cp$dat)
  ggsave(file.path(fig_dir, "FigS6_Reverse_CogPerf_Scatter.pdf"),
         p_rev2[[1]] + theme_nature + ggtitle("Reverse MR: HF on Cognitive Performance"), width = 7, height = 6.5)
}

cat("\n所有图表已保存至", fig_dir, "\n")
