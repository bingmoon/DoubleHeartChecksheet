# ============================================================
# E01_import_data.R
# 下载并清洗扩展训练集：NHANES 2009-2018 (F, G, H, I, J)
# 输出：extended_cohort/Extended_Cohort_Data.RData
# ============================================================

library(nhanesA)
library(dplyr)
library(tidyr)
library(readr)

ext_dir <- "extended_cohort"
if (!dir.exists(ext_dir)) dir.create(ext_dir, recursive = TRUE)

# 解析函数
parse_mcq <- function(x) {
  ifelse(x == 1 | grepl("Yes", as.character(x), ignore.case = TRUE), 1, 0)
}

parse_dpq <- function(x) {
  char_x <- as.character(x)
  case_when(
    char_x %in% c("0") | grepl("Not at all", char_x, ignore.case = TRUE) ~ 0,
    char_x %in% c("1") | grepl("Several days", char_x, ignore.case = TRUE) ~ 1,
    char_x %in% c("2") | grepl("More than half", char_x, ignore.case = TRUE) ~ 2,
    char_x %in% c("3") | grepl("Nearly every day", char_x, ignore.case = TRUE) ~ 3,
    TRUE ~ NA_real_
  )
}

# 五个周期：2009-2018 全覆盖
cycles <- c("F", "G", "H", "I", "J")
all_raw <- list()

for (cyc in cycles) {
  cat("Downloading cycle", cyc, "...\n")
  
  demo <- nhanes(paste0("DEMO_", cyc)) %>%
    select(SEQN, RIAGENDR, RIDAGEYR)
  dpq  <- nhanes(paste0("DPQ_", cyc)) %>%
    select(SEQN, DPQ010, DPQ020, DPQ030, DPQ040, DPQ050, DPQ060, DPQ070, DPQ080, DPQ090)
  mcq  <- nhanes(paste0("MCQ_", cyc)) %>%
    select(SEQN, MCQ160B, MCQ160C, MCQ160E, MCQ160F)
  bmx  <- nhanes(paste0("BMX_", cyc)) %>%
    select(SEQN, BMXBMI)
  
  raw <- demo %>%
    left_join(dpq, by = "SEQN") %>%
    left_join(mcq, by = "SEQN") %>%
    left_join(bmx, by = "SEQN")
  
  mort_file <- switch(cyc,
    "F" = "NHANES_2009_2010_MORT_2019_PUBLIC.dat",
    "G" = "NHANES_2011_2012_MORT_2019_PUBLIC.dat",
    "H" = "NHANES_2013_2014_MORT_2019_PUBLIC.dat",
    "I" = "NHANES_2015_2016_MORT_2019_PUBLIC.dat",
    "J" = "NHANES_2017_2018_MORT_2019_PUBLIC.dat"
  )
  if (!file.exists(mort_file)) {
    download.file(paste0("https://ftp.cdc.gov/pub/health_statistics/NCHS/datalinkage/linked_mortality/", mort_file),
                  destfile = mort_file, mode = "wb")
  }
  
  ndi_positions <- fwf_positions(
    start = c(1, 15, 16, 17, 43),
    end   = c(14, 15, 16, 19, 45),
    col_names = c("SEQN", "ELIGSTAT", "MORTSTAT", "UCOD_LEADING", "PERMTH_INT")
  )
  
  mort <- read_fwf(mort_file, col_positions = ndi_positions, show_col_types = FALSE) %>%
    mutate(SEQN = as.numeric(SEQN),
           PERMTH_INT = as.numeric(PERMTH_INT),
           MORTSTAT = as.numeric(MORTSTAT))
  
  raw <- raw %>% left_join(mort, by = "SEQN")
  all_raw[[cyc]] <- raw
}

# 合并与清洗
DH_Data_Ext <- bind_rows(all_raw) %>%
  mutate(
    across(c(MCQ160B, MCQ160C, MCQ160E, MCQ160F), parse_mcq),
    Has_CVD = ifelse(MCQ160B == 1 | MCQ160C == 1 | MCQ160E == 1 | MCQ160F == 1, 1, 0)
  ) %>%
  filter(Has_CVD == 1, !is.na(MORTSTAT), !is.na(PERMTH_INT)) %>%
  mutate(across(starts_with("DPQ"), parse_dpq)) %>%
  drop_na(starts_with("DPQ")) %>%
  mutate(
    Cognitive_Score = DPQ010 + DPQ020 + DPQ060 + DPQ070 + DPQ090,
    PHQ9_Total = rowSums(select(., starts_with("DPQ"))),
    Gender = factor(ifelse(as.character(RIAGENDR) %in% c("2", "Female"), "Female", "Male")),
    Age = as.numeric(RIDAGEYR),
    BMI = as.numeric(BMXBMI),
    Survival_Months = PERMTH_INT,
    CVD_Death = ifelse(MORTSTAT == 1 & grepl("001|002|003|004", as.character(UCOD_LEADING)), 1, 0)
  ) %>%
  filter(Age >= 40, !is.na(BMI), Survival_Months > 0)

save(DH_Data_Ext, file = file.path(ext_dir, "Extended_Cohort_Data.RData"))

cat(sprintf("\n扩展队列样本量: %d, 事件数: %d\n", nrow(DH_Data_Ext), sum(DH_Data_Ext$CVD_Death)))
