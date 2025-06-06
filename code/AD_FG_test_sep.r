# Databricks notebook source
library(foreach)
install.packages("doParallel")
library(doParallel)
install.packages("cmprsk")
library(cmprsk)
library(tidyverse)
library(lubridate)
raw_data_path = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/mingzhoufu@mednet.ucla.edu/AD_trajectory/data/"
# Register the parallel backend and specify the number of cores
nCores = detectCores() - 1  # Use one less than the total number of cores
print(nCores)
registerDoParallel(cores = nCores)

# COMMAND ----------

load(file = paste0(raw_data_path, "FG_test/AD_risk/combined_visit_wide_short_0610.rda"))
load(file = paste0(raw_data_path, "upload_ucla/embed_sim_cal.rda"))
load(file = paste0(raw_data_path, "FG_test/AD_risk/FG_AD_sample_demo_0610.rda"))
embed_icd = colnames(embed_sim_cal)
print(length(embed_icd)) # 1619
test_icd = intersect(names(combined_visit_wide_short), embed_icd) %>% sort()
print(length(test_icd)) # 462
# also remove certain chapter codes
rm_codes = grep("^[OPQRVXY][0-9]+", test_icd, value = TRUE)
icd_lst_filter = setdiff(test_icd, rm_codes) 
print(length(icd_lst_filter)) # N = 388

# COMMAND ----------

table(FG_AD_sample_demo$outcome_status)

# COMMAND ----------

FG_result = foreach(i = 1:50, .combine = 'rbind', .packages = c("dplyr", "survival")) %dopar% {
  print(i)
  exposure_icd = icd_lst_filter[i]
  exposure_label = combined_visit_wide_short %>% select(PatientID, all_of(exposure_icd)) %>% rename("exposure_status" = exposure_icd)
  full_test_df = FG_AD_sample_demo %>% inner_join(exposure_label)
  # check for EPV
  event_per_variable = table(full_test_df$exposure_status, full_test_df$outcome_status)
  if (all(event_per_variable >= 10)) {
    # prepare dataset
    cov_matrix = model.matrix(~ factor(exposure_status) + age_first_visit + female + factor(race_ethnicity), data = full_test_df)[, -1]
    crr_model = crr(ftime = full_test_df$follow_up_time, fstatus = full_test_df$outcome_status, cov1 = cov_matrix, failcode = 1, cencode = 0)
    summary_crr = summary(crr_model)
    coefficient = summary_crr$coef[1,] %>% t() %>% as.data.frame() %>% mutate(exposure_ICD = exposure_icd, outcome_ICD = "G30")
  } else {
    coefficient = c(rep(NA, 5), exposure_icd, outcome_icd) %>% as.data.frame() %>% t()
    colnames(coefficient) = c("coef", "exp(coef)", "se(coef)", "z", "p-value", "exposure_ICD", "outcome_ICD")
  }
  coefficient  # The last expression in the loop is returned and combined
}

# COMMAND ----------

save(FG_result, file = paste0(raw_data_path, "FG_test/AD_risk/FG_result/FG_result_1.rda"))

# COMMAND ----------

head(FG_result)