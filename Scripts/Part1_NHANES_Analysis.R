# ==============================================================================
# 项目名称: Double Heart (DH)
# 分析模块: 第一部分：NHANES 真实世界数据挖掘与生存分析
# 步骤名称: Step 1: NHANES 核心数据流拉取、多周期合并与死亡随访数据匹配 (修正版)
# ==============================================================================

# 1. 环境初始化与目录设置 
setwd("/Users/bing/DH")

# 检查并自动安装必备包
req_packages <- c("nhanesA", "dplyr", "tidyr", "readr")
new_packages <- req_packages[!(req_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(nhanesA)
library(dplyr)
library(tidyr)
library(readr)

# 延长 R 网络超时时间
options(timeout = 600)

# ==============================================================================
# 2. 从 NHANES 数据库拉取基线问卷与检验数据
# ==============================================================================

fetch_nhanes_cycle <- function(cycle_letter) {
  cat(paste0(">>> 正在通过 nhanesA 接口抓取周期: ", cycle_letter, " 的基线数据...\n"))
  
  # A. 人口学数据
  demo <- nhanes(paste0("DEMO_", cycle_letter)) %>% 
    select(SEQN, RIAGENDR, RIDAGEYR, DMDEDUC2, INDFMPIR) 
  
  # B. 抑郁症筛查问卷 (PHQ-9)
  dpq <- nhanes(paste0("DPQ_", cycle_letter)) %>%
    select(SEQN, DPQ010, DPQ020, DPQ030, DPQ040, DPQ050, DPQ060, DPQ070, DPQ080, DPQ090)
  
  # C. 医疗状况 (心血管病史)
  mcq <- nhanes(paste0("MCQ_", cycle_letter)) %>%
    select(SEQN, MCQ160B, MCQ160C, MCQ160E, MCQ160F) 
  
  # D. 身体测量 (BMI)
  bmx <- nhanes(paste0("BMX_", cycle_letter)) %>%
    select(SEQN, BMXBMI)
  
  merged_cycle <- demo %>%
    left_join(dpq, by = "SEQN") %>%
    left_join(mcq, by = "SEQN") %>%
    left_join(bmx, by = "SEQN")
  
  return(merged_cycle)
}

# 合并 2015-2016 (I) 和 2017-2018 (J) 
cohort_I <- fetch_nhanes_cycle("I")
cohort_J <- fetch_nhanes_cycle("J")
raw_cohort <- bind_rows(cohort_I, cohort_J)

cat(">>> 基线数据拉取完毕，开始获取生存随访数据...\n")

# ==============================================================================
# 3. 解析 CDC NDI 固定宽度 (.dat) 死亡登记数据 (修复数据类型冲突)
# ==============================================================================

# 定义 CDC 官方 .dat 文件的列宽 (固定宽度截取)
mort_fwf_cols <- fwf_widths(
  c(14, 1, 1, 3, 1, 1, 3, 3), 
  c("SEQN", "ELIGSTAT", "MORTSTAT", "UCOD_LEADING", "DIABETES", "HYPERTEN", "PERMTH_INT", "PERMTH_EXM")
)

# 读入死亡数据，并强制统一 SEQN 为数值型 (as.numeric)
cat(">>> 正在解析死亡登记 .dat 数据并统一数据类型...\n")
mort_2015 <- read_fwf(file_2015, col_positions = mort_fwf_cols, show_col_types = FALSE) %>%
  mutate(SEQN = as.numeric(SEQN))

mort_2017 <- read_fwf(file_2017, col_positions = mort_fwf_cols, show_col_types = FALSE) %>%
  mutate(SEQN = as.numeric(SEQN))

# 现在可以安全合并了
mortality_data <- bind_rows(mort_2015, mort_2017) %>%
  select(SEQN, MORTSTAT, PERMTH_INT, UCOD_LEADING)

# ==============================================================================
# 4. 构建最终分析大表并保存
# ==============================================================================

cat(">>> 正在进行基线特征与生存结局的内连接 (Inner Join)...\n")
DH_Step1_Data <- raw_cohort %>%
  inner_join(mortality_data, by = "SEQN")

saveRDS(DH_Step1_Data, file = "DH_Step1_RawData.rds")

cat("\n==============================================================================\n")
cat("SUCCESS! 第一部分 Step 1 终于圆满执行完毕！\n")
cat(paste0("总样本量: ", nrow(DH_Step1_Data), " 例\n"))
cat("数据已安全打包保存为: /Users/bing/DH/DH_Step1_RawData.rds\n")
cat("==============================================================================\n")

# ==============================================================================
# 项目名称: Double Heart (DH)
# 分析模块: 第一部分：NHANES 真实世界数据挖掘与生存分析
# 步骤名称: Step 2: 目标心血管队列筛选与特征重构 (全容错工业级版)
# ==============================================================================

# 1. 环境初始化与目录设置
setwd("/Users/bing/DH")
library(dplyr)
library(tidyr)
library(readr)

cat(">>> 1. 正在使用【绝对坐标系】精准提取 CDC 死亡登记数据...\n")
# 抛弃容易出错的列宽相加，直接指定 2019 NDI 字典的绝对起始和终止位置
# 完美跨过被抹除的 22-42 列隐私数据区
ndi_positions <- fwf_positions(
  start = c(1, 15, 16, 17, 43),
  end   = c(14, 15, 16, 19, 45),
  col_names = c("SEQN", "ELIGSTAT", "MORTSTAT", "UCOD_LEADING", "PERMTH_INT")
)

# 极速读取本地已有的 .dat 文件
mort_2015_safe <- read_fwf("NHANES_2015_2016_MORT_2019_PUBLIC.dat", col_positions = ndi_positions, show_col_types = FALSE) %>% mutate(SEQN = as.numeric(SEQN))
mort_2017_safe <- read_fwf("NHANES_2017_2018_MORT_2019_PUBLIC.dat", col_positions = ndi_positions, show_col_types = FALSE) %>% mutate(SEQN = as.numeric(SEQN))

mortality_safe <- bind_rows(mort_2015_safe, mort_2017_safe) %>%
  mutate(
    PERMTH_INT = as.numeric(PERMTH_INT),
    MORTSTAT = as.numeric(MORTSTAT)
  )

cat(">>> 2. 正在加载基线库并无缝融合准确的生存数据...\n")
# 读取 Step 1 存好的大表，清洗掉旧的死亡列，挂载新提取的安全列
DH_Step1_Data <- readRDS("DH_Step1_RawData.rds") %>%
  select(-any_of(c("MORTSTAT", "PERMTH_INT", "UCOD_LEADING", "ELIGSTAT"))) %>%
  left_join(mortality_safe, by = "SEQN")


# ==============================================================================
# 构建智能解析器：击碎标签与缺失值陷阱
# ==============================================================================
# 解析心血管病史 (兼容数字1与文本"Yes")
parse_mcq <- function(x) {
  ifelse(x == 1 | grepl("Yes", as.character(x), ignore.case = TRUE), 1, 0)
}

# 解析PHQ-9抑郁条目 (兼容0-3数字与NHANES英文问卷标签，自动转译)
parse_dpq <- function(x) {
  char_x <- as.character(x)
  case_when(
    char_x %in% c("0") | grepl("Not at all", char_x, ignore.case = TRUE) ~ 0,
    char_x %in% c("1") | grepl("Several days", char_x, ignore.case = TRUE) ~ 1,
    char_x %in% c("2") | grepl("More than half", char_x, ignore.case = TRUE) ~ 2,
    char_x %in% c("3") | grepl("Nearly every day", char_x, ignore.case = TRUE) ~ 3,
    TRUE ~ NA_real_ # 将拒答、不知道等全部强制转为 NA
  )
}

cat(">>> 3. 正在执行心血管病变人群精准筛选与抑郁维度重构...\n")

DH_Step2_Data <- DH_Step1_Data %>%
  # A. 解析合并症并提取 CVD 队列
  mutate(
    CVD_B = parse_mcq(MCQ160B),
    CVD_C = parse_mcq(MCQ160C),
    CVD_E = parse_mcq(MCQ160E),
    CVD_F = parse_mcq(MCQ160F),
    Has_CVD = ifelse(CVD_B == 1 | CVD_C == 1 | CVD_E == 1 | CVD_F == 1, 1, 0)
  ) %>%
  filter(Has_CVD == 1) %>% 
  
  # B. 剔除死亡状态(MORTSTAT)或生存月数(PERMTH_INT)不明的样本
  filter(!is.na(MORTSTAT) & !is.na(PERMTH_INT)) %>%
  
  # C. 标准化 PHQ-9 所有条目并剔除问卷不完整的样本
  mutate(across(starts_with("DPQ"), parse_dpq)) %>%
  drop_na(starts_with("DPQ")) %>%
  
  # D. 核心变量重构：双维度得分计算
  mutate(
    Somatic_Score = DPQ030 + DPQ040 + DPQ050 + DPQ080,
    Cognitive_Score = DPQ010 + DPQ020 + DPQ060 + DPQ070 + DPQ090,
    PHQ9_Total = Somatic_Score + Cognitive_Score
  ) %>%
  
  # E. 整理最终模型字段
  mutate(
    Gender = factor(ifelse(RIAGENDR == 1 | grepl("Male", as.character(RIAGENDR), ignore.case = TRUE), "Male", "Female")),
    Age = as.numeric(RIDAGEYR),
    BMI = as.numeric(BMXBMI),
    Survival_Months = PERMTH_INT,
    AllCause_Death = MORTSTAT,
    # 提取心因性死亡: CDC UCOD_LEADING 含有 001,002,003,004 (心脏相关) 时记为 1
    CVD_Death = ifelse(MORTSTAT == 1 & grepl("001|002|003|004", as.character(UCOD_LEADING)), 1, 0)
  ) %>%
  
  # 剔除 年龄<40岁 及 缺失 BMI 的样本
  filter(Age >= 40 & !is.na(BMI)) %>%
  
  # F. 数据瘦身
  select(SEQN, Gender, Age, BMI, starts_with("DPQ"), 
         Somatic_Score, Cognitive_Score, PHQ9_Total, 
         Survival_Months, AllCause_Death, CVD_Death)

saveRDS(DH_Step2_Data, file = "DH_Step2_CleanedData.rds")

cat("\n==============================================================================\n")
cat("SUCCESS! 第一部分 Step 2 (全容错版) 执行完毕！\n")
cat(paste0("有效 CVD 队列样本量: ", nrow(DH_Step2_Data), " 例\n"))
cat("死亡事件数 (全因): ", sum(DH_Step2_Data$AllCause_Death == 1, na.rm=TRUE), " 例\n")
cat("心血管专病死亡数: ", sum(DH_Step2_Data$CVD_Death == 1, na.rm=TRUE), " 例\n")
cat("==============================================================================\n")

# ==============================================================================
# 项目名称: Double Heart (DH)
# 分析模块: 第一部分：NHANES 真实世界数据挖掘与生存分析
# 步骤名称: Step 3: 单因素 Cox 预后效能排名与 Top-K 核心症状筛选 (强效挽救版)
# ==============================================================================

# 1. 环境初始化
setwd("/Users/bing/DH")

req_packages <- c("dplyr", "survival", "ggplot2")
new_packages <- req_packages[!(req_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(dplyr)
library(survival)
library(ggplot2)

cat(">>> 1. 正在加载心血管确诊队列数据...\n")
DH_Data <- readRDS("DH_Step2_CleanedData.rds")

# ==============================================================================
# 2. 遍历执行单因素 Cox 比例风险回归
# ==============================================================================
cat(">>> 2. 正在让 9 个抑郁条目分别与心血管死亡进行单因素 Cox 决斗...\n")

# 定义 9 个 PHQ-9 条目
dpq_vars <- c("DPQ010", "DPQ020", "DPQ030", "DPQ040", "DPQ050", "DPQ060", "DPQ070", "DPQ080", "DPQ090")

# 建立空数据框储存结果
ranking_results <- data.frame()

for (var in dpq_vars) {
  # 构建单因素生存公式
  formula_str <- paste("Surv(Survival_Months, CVD_Death) ~", var)
  fit <- coxph(as.formula(formula_str), data = DH_Data)
  
  # 提取统计量
  s <- summary(fit)
  
  # 如果该条目模型成功拟合，记录其危险系数和 Z-score
  if (nrow(s$coefficients) > 0) {
    ranking_results <- rbind(ranking_results, data.frame(
      Symptom = var,
      HR = s$coefficients[1, "exp(coef)"],
      Z_Score = abs(s$coefficients[1, "z"]),  # Z-score的绝对值代表预测效力的大小
      P_Value = s$coefficients[1, "Pr(>|z|)"]
    ))
  }
}

# ==============================================================================
# 3. 计算预测排名与 Top-K 特征提取
# ==============================================================================
cat(">>> 3. 正在根据 Z-score (预测致命性) 进行排名，并锁定 Top 3 靶点...\n")

# 按 Z-Score 从大到小排序
ranking_results <- ranking_results[order(-ranking_results$Z_Score), ]

# 添加临床易读的症状标签 (为作图准备)
symptom_labels <- c(
  "DPQ010" = "1. Anhedonia (兴趣缺失)", "DPQ020" = "2. Depressed Mood (抑郁情绪)",
  "DPQ030" = "3. Sleep Disturbance (睡眠障碍)", "DPQ040" = "4. Fatigue (疲劳)",
  "DPQ050" = "5. Appetite Change (食欲改变)", "DPQ060" = "6. Guilt/Worthlessness (自责)",
  "DPQ070" = "7. Concentration Issues (注意力差)", "DPQ080" = "8. Psychomotor Agitation (精神运动异常)",
  "DPQ090" = "9. Suicidal Thoughts (自杀意念)"
)
ranking_results$Label <- symptom_labels[ranking_results$Symptom]

# 强制提取排名前 3 的症状作为后续多因素模型的核心特征
top_k <- 3
selected_features <- ranking_results$Symptom[1:top_k]

cat("========================================================\n")
cat("🔥 挽救成功！最具心血管致死预测效力的 Top 3 症状为:\n")
for(i in 1:top_k) {
  cat(sprintf("第 %d 名: %s (Z-score: %.2f, P-value: %.3f)\n", 
              i, ranking_results$Label[i], ranking_results$Z_Score[i], ranking_results$P_Value[i]))
}
cat("========================================================\n")

# 将选中的特征覆盖保存至 RDS，供 Step 4 提取
saveRDS(selected_features, file = "DH_Step3_SelectedFeatures.rds")

# ==============================================================================
# 4. 生成 SCI 级别高精度 PDF：预后效能棒棒糖图 (Lollipop Chart)
# ==============================================================================
cat(">>> 4. 正在渲染并导出 SCI 级别高精度特征重要性棒棒糖图...\n")

# 将数据框排序锁定为因子，保证绘图顺序
ranking_results$Label <- factor(ranking_results$Label, levels = rev(ranking_results$Label))

# 区分被选中的 Top 3 和落选条目 (赋予不同颜色)
ranking_results$Selected <- ifelse(ranking_results$Symptom %in% selected_features, "Selected (Top 3)", "Unselected")

# 绘制极具设计感的棒棒糖图
lollipop_plot <- ggplot(ranking_results, aes(x = Label, y = Z_Score, color = Selected)) +
  geom_segment(aes(x = Label, xend = Label, y = 0, yend = Z_Score), size = 1.2) +
  geom_point(size = 5) +
  scale_color_manual(values = c("Selected (Top 3)" = "#D62728", "Unselected" = "#A9A9A9")) +
  coord_flip() + # 翻转坐标轴，方便阅读长文本
  theme_classic() +
  labs(title = "Predictive Power of PHQ-9 Symptoms for Cardiovascular Mortality",
       subtitle = "Ranked by Univariable Cox Wald Z-scores",
       x = "PHQ-9 Depression Items",
       y = "Predictive Power (Absolute Z-score)") +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    axis.text.x = element_text(size = 11),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# 使用 ggsave 安全导出
ggsave(filename = "DH_Part1_Step3_FeatureRanking_Lollipop.pdf", 
       plot = lollipop_plot, width = 8.5, height = 6, device = "pdf")

cat("SUCCESS! 第一部分 Step 3 (强效挽救版) 执行完毕！\n")
cat("已生成极其优雅的特征重要性排名图: DH_Part1_Step3_FeatureRanking_Lollipop.pdf\n")

# ==============================================================================
# 项目名称: Double Heart (DH)
# 分析模块: 第一部分：NHANES 真实世界数据挖掘与生存分析
# 步骤名称: Step 4: 多因素 Cox 模型构建、结果导出与高级可视化 (终极全功能版)
# ==============================================================================

# 1. 环境初始化与目录设置
setwd("/Users/bing/DH")

# 检查并自动安装必备包
req_packages <- c("dplyr", "survival", "survminer", "ggplot2")
new_packages <- req_packages[!(req_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(dplyr)
library(survival)
library(survminer)
library(ggplot2)

cat(">>> 1. 正在载入队列数据与 Top 3 核心心理特征...\n")
Raw_Data <- readRDS("DH_Step2_CleanedData.rds")
selected_features <- readRDS("DH_Step3_SelectedFeatures.rds")

# ==============================================================================
# 2. 自动化数据校验：内嵌性别分类精确修正逻辑
# ==============================================================================
cat(">>> 2. 正在进行基线特征数据清洗与性别分层校验...\n")
Step1_Raw <- readRDS("DH_Step1_RawData.rds") %>% select(SEQN, RIAGENDR)

DH_Data <- Raw_Data %>%
  select(-any_of("Gender")) %>% 
  left_join(Step1_Raw, by = "SEQN") %>%
  mutate(
    # 严格根据原始编码区分女性，彻底斩断字符包含引起的分类错误
    Gender = factor(ifelse(as.character(RIAGENDR) %in% c("2", "Female"), "Female", "Male"),
                    levels = c("Female", "Male"))
  ) %>%
  select(-RIAGENDR)

# ==============================================================================
# 3. 构建多因素 Cox 比例风险回归模型
# ==============================================================================
cat(">>> 3. 正在构建校正了年龄、性别、BMI 的多因素 Cox 生存回归模型...\n")

# 动态组装高级生存回归公式
cox_formula_str <- paste("Surv(Survival_Months, CVD_Death) ~", 
                         paste(selected_features, collapse = " + "), 
                         "+ Age + Gender + BMI")
cox_formula <- as.formula(cox_formula_str)

# 拟合 Cox 模型
cox_model <- coxph(cox_formula, data = DH_Data)

cat("\n=================== 多因素 Cox 模型回归摘要 ===================\n")
print(summary(cox_model))
cat("===============================================================\n")

# ==============================================================================
# 4. [新增] 提取 Cox 模型核心结果并一键导出为 CSV 论文表格
# ==============================================================================
cat(">>> 4. 正在提取 Cox 模型的 HR 值、95% CI 和 P 值，并导出为 CSV...\n")

cox_summary <- summary(cox_model)

# 组装高标准的数据框
cox_results <- data.frame(
  Variable = rownames(cox_summary$coefficients),
  Hazard_Ratio = round(cox_summary$coefficients[, "exp(coef)"], 3),
  CI_Lower_95 = round(cox_summary$conf.int[, "lower .95"], 3),
  CI_Upper_95 = round(cox_summary$conf.int[, "upper .95"], 3),
  P_Value = signif(cox_summary$coefficients[, "Pr(>|z|)"], 3)
)

# 增加一列直观的显著性星号标记
cox_results$Significance <- ifelse(cox_results$P_Value < 0.001, "***",
                                   ifelse(cox_results$P_Value < 0.01, "**",
                                          ifelse(cox_results$P_Value < 0.05, "*", "")))

# 输出保存为 CSV，供后续 Word/Excel 直接排版
write.csv(cox_results, file = "DH_Part1_Step4_Cox_Results.csv", row.names = FALSE)

# ==============================================================================
# 5. 绘制并导出森林图 (使用 ggsave 彻底消灭第一页空白 Bug)
# ==============================================================================
cat(">>> 5. 正在通过 ggsave 渲染并导出高分辨率 Cox 危险系数森林图...\n")

forest_plot <- ggforest(cox_model, data = DH_Data, 
                        main = "Hazard Ratios for Cardiovascular Mortality",
                        cpositions = c(0.02, 0.22, 0.4),
                        fontsize = 1.0, 
                        refLabel = "reference", noDigits = 2)

ggsave(filename = "DH_Part1_Step4_Cox_ForestPlot.pdf", 
       plot = forest_plot, 
       width = 8, 
       height = 6, 
       device = "pdf")

# ==============================================================================
# 6. 基于模型预测值建立综合心理风险预警评分 (Risk Score)
# ==============================================================================
cat(">>> 6. 正在计算个体生存风险指数并执行中位数黄金切分...\n")

DH_Data$Risk_Score <- predict(cox_model, type = "risk")

median_risk <- median(DH_Data$Risk_Score, na.rm = TRUE)
DH_Data$Risk_Group <- factor(ifelse(DH_Data$Risk_Score > median_risk, "High Risk", "Low Risk"),
                             levels = c("Low Risk", "High Risk"))

# ==============================================================================
# 7. 绘制并导出 Kaplan-Meier 生存曲线
# ==============================================================================
cat(">>> 7. 正在渲染并导出 Kaplan-Meier 累积生存概率曲线...\n")

km_fit <- survfit(Surv(Survival_Months, CVD_Death) ~ Risk_Group, data = DH_Data)

pdf("DH_Part1_Step4_KM_SurvivalCurve.pdf", width = 7.5, height = 6.5, onefile = FALSE)

km_plot <- ggsurvplot(km_fit,
                      data = DH_Data,
                      pval = TRUE,             
                      pval.coord = c(10, 0.15),
                      conf.int = TRUE,         
                      risk.table = TRUE,       
                      risk.table.col = "strata", 
                      risk.table.height = 0.24,
                      palette = c("#1F77B4", "#D62728"), # 科技蓝与医学红经典配色
                      title = "Kaplan-Meier Curve for Cardiovascular Survival",
                      xlab = "Follow-up Time (Months)", 
                      ylab = "Survival Probability",
                      legend.title = "Risk Stratification",
                      legend.labs = c("Low Psychological Risk", "High Psychological Risk"),
                      ggtheme = theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold")))

print(km_plot)
dev.off()

# ==============================================================================
# 8. 最终分析大表持久化备份
# ==============================================================================
saveRDS(DH_Data, file = "DH_Step4_FinalData_with_Risk.rds")

cat("\n==============================================================================\n")
cat("SUCCESS! 第一部分 Step 4 (含 CSV 导出) 震撼执行完毕！\n")
cat("请在 /Users/bing/DH 目录下查收您的科研成果:\n")
cat(" 1. DH_Part1_Step4_Cox_Results.csv (核心数据表格，可直接导入 Excel)\n")
cat(" 2. DH_Part1_Step4_Cox_ForestPlot.pdf (单页完美森林图)\n")
cat(" 3. DH_Part1_Step4_KM_SurvivalCurve.pdf (高辨识度双色生存曲线)\n")
cat("==============================================================================\n")

# ==============================================================================
# 项目名称: Double Heart (DH)
# 补充模块: 生成 SCI 论文级基线特征表 (Table 1)
# ==============================================================================

# 1. 环境初始化
setwd("/Users/bing/DH")

# 检查并安装强大的 tableone 包（专为医学论文 Table 1 设计）
if (!require("tableone")) install.packages("tableone")
library(tableone)
library(dplyr)

cat(">>> 1. 正在载入带有风险分层的终极队列数据...\n")
DH_Data <- readRDS("DH_Step4_FinalData_with_Risk.rds")

# 2. 定义 Table 1 需要展示的变量
# 包含：人口学指标、代谢指标、核心心理预警靶点得分，以及生存结局
vars_to_summarize <- c("Age", "Gender", "BMI", 
                       "DPQ070", "DPQ030", "DPQ040", 
                       "PHQ9_Total", "Survival_Months", "CVD_Death")

# 定义哪些变量是分类变量 (Categorical variables)
categorical_vars <- c("Gender", "CVD_Death")

# 3. 构建 Table 1 对象
# 这里我们选择按照 "Risk_Group" (低风险 vs 高风险) 进行分组对比，这是顶刊最喜欢的展示方式
cat(">>> 2. 正在按高低心理风险分组计算基线特征与 P 值...\n")
table1_obj <- CreateTableOne(vars = vars_to_summarize, 
                             factorVars = categorical_vars, 
                             strata = "Risk_Group", 
                             data = DH_Data, 
                             addOverall = TRUE) # 同时保留整体队列 (Overall) 的统计

# 4. 格式化输出
# 连续变量非正态分布时，输出中位数和四分位距 (非必须，这里按标准均值±标准差输出)
table1_printed <- print(table1_obj, 
                        showAllLevels = TRUE, 
                        quote = FALSE, 
                        noSpaces = TRUE, 
                        printToggle = FALSE)

# 5. 导出为 CSV 文件供直接插入论文
write.csv(table1_printed, file = "DH_Part1_Table1_BaselineCharacteristics.csv")

cat("==============================================================================\n")
cat("SUCCESS! Table 1 (基线特征表) 已成功生成！\n")
cat("文件已保存为: /Users/bing/DH/DH_Part1_Table1_BaselineCharacteristics.csv\n")
cat("您可以直接使用 Excel 打开并排版到您的 LaTeX 文档中。\n")
cat("==============================================================================\n")
