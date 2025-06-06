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
'%!in%' = function(x,y)!{'%in%'(x,y)}

# COMMAND ----------

load(file = paste0(raw_data_path, "FG_test/other_risk/diag_beforeAD_filter_0610.rda"))
load(file = paste0(raw_data_path, "patient_data/mod/AD_demographic_0610.rda"))
print(dim(AD_demographic)) # dim = (24473,15)
load(file = paste0(raw_data_path, "patient_data/mod/age_first_visit_add_0610.rda"))
print(dim(age_first_visit_add)) # dim = (919590,2)
load(file = paste0(raw_data_path, "upload_ucla/embed_sim_cal.rda"))
embed_icd = colnames(embed_sim_cal)
print(length(embed_icd)) # 1619
load(file = paste0(raw_data_path, "FG_test/other_risk/icd_lst_final_0610.rda"))
print(length(icd_lst_final)) # 302

# COMMAND ----------

outcome_icd = "F02"

# COMMAND ----------

# identify patients with outcome ICD
outcome_define = diag_beforeAD_filter %>% filter(ICD_3digit == outcome_icd) %>%
  select(PatientID, EncounterDate) %>% rename(OutcomeDate = EncounterDate) %>% unique()
outcome_patenc_filter = diag_beforeAD_filter %>% filter(PatientID %in% outcome_define$PatientID) %>% 
  left_join(outcome_define) %>% filter(EncounterDate <= OutcomeDate)
# define outcome status
status_label = AD_demographic %>% 
  mutate(outcome_status = case_when(
    PatientID %!in% outcome_define$PatientID & PatientLivingStatus == "Deceased" ~ 2,
    PatientID %in% outcome_define$PatientID ~ 1,
    TRUE ~ 0
  )) %>% select(PatientID, outcome_status)
print(table(status_label$outcome_status))
# calculate new EHR features
follow_up_time_case = outcome_patenc_filter %>% group_by(PatientID) %>%
  summarise(follow_up_time = as.numeric(OutcomeDate - min(EncounterDate))/365.25) %>% ungroup() %>% 
  select(PatientID, follow_up_time) %>% unique()
follow_up_time_control = AD_demographic %>% filter(PatientID %!in% outcome_define$PatientID) %>%
  select(PatientID, record_length) %>% rename("follow_up_time" = "record_length")
follow_up_time_df = rbind(follow_up_time_case, follow_up_time_control)
# finalize sample demographic covariates
sample_demo = AD_demographic %>% left_join(age_first_visit_add) %>% left_join(status_label) %>%
  left_join(follow_up_time_df) %>% select(PatientID, age_first_visit, female, race_ethnicity, follow_up_time, outcome_status) %>% unique()
case_visit_wide = outcome_patenc_filter %>%
  select(PatientID, ICD_3digit) %>% mutate(status = 1) %>% unique() %>%
  pivot_wider(names_from = ICD_3digit, values_from = status, values_fill = 0) %>% 
  mutate_all(~ifelse(is.na(.), 0, .)) 
# combine all together to make a big dummy dataset
control_visit_wide = diag_beforeAD_filter %>% filter(PatientID %!in% outcome_define$PatientID) %>%
  select(PatientID, ICD_3digit) %>% mutate(status = 1) %>% unique() %>%
  pivot_wider(names_from = ICD_3digit, values_from = status, values_fill = 0) %>% 
  mutate_all(~ifelse(is.na(.), 0, .)) 
# pull ICD list
outcome_step_to_test = outcome_patenc_filter %>% pull(ICD_3digit) %>% unique() %>% sort()
icd_to_test = setdiff(outcome_step_to_test, outcome_icd)
icd_to_test = intersect(icd_to_test, names(control_visit_wide))
print(paste0("Length of ICD to test: ", length(icd_to_test))) # 300
control_visit_wide_short = control_visit_wide %>% select("PatientID", all_of(icd_to_test))
case_visit_wide_short = case_visit_wide %>% select("PatientID", all_of(icd_to_test))
combined_visit_wide = rbind(case_visit_wide_short, control_visit_wide_short)
print(dim(combined_visit_wide)) # dim = (24473,301)

# COMMAND ----------

# Example names of the objects to keep
keep_var = c("icd_to_test", "combined_visit_wide", "nCores", "sample_demo", "outcome_icd", "raw_data_path")
# List all objects in the environment
all_objects = ls()
# Remove all data frames except the one specified
for (obj_name in all_objects) {
  if (obj_name %!in% keep_var) {
    if (is.data.frame(get(obj_name))) {
      rm(list = obj_name, envir = .GlobalEnv)
    }
  }
}
# Trigger garbage collection
gc()

# COMMAND ----------

# Use foreach for the parallel loop
FG_result = foreach(i = 1:length(icd_to_test), .combine = 'rbind', .packages = c("dplyr", "survival")) %dopar% {
  print(i)
  exposure_icd = icd_to_test[i]
  exposure_label = combined_visit_wide %>% select(PatientID, all_of(exposure_icd)) %>% rename("exposure_status" = exposure_icd)
  full_test_df = sample_demo %>% inner_join(exposure_label)
  # check for EPV
  event_per_variable = table(full_test_df$exposure_status, full_test_df$outcome_status)
  if (all(event_per_variable >= 10)) {
    # prepare dataset
    cov_matrix = model.matrix(~ factor(exposure_status) + age_first_visit + female + factor(race_ethnicity), data = full_test_df)[, -1]
    crr_model = crr(ftime = full_test_df$follow_up_time, fstatus = full_test_df$outcome_status, cov1 = cov_matrix, failcode = 1, cencode = 0)
    summary_crr = summary(crr_model)
    coefficient = summary_crr$coef[1,] %>% t() %>% as.data.frame() %>% mutate(exposure_ICD = exposure_icd, outcome_ICD = outcome_icd)
  } else {
    coefficient = c(rep(NA, 5), exposure_icd, outcome_icd) %>% as.data.frame() %>% t()
    colnames(coefficient) = c("coef", "exp(coef)", "se(coef)", "z", "p-value", "exposure_ICD", "outcome_ICD")
  }
  coefficient  # The last expression in the loop is returned and combined
  save(coefficient, file = paste0(raw_data_path, "FG_test/other_risk/FG_result/", outcome_icd, "/", exposure_icd, ".rda"))
}

# COMMAND ----------

