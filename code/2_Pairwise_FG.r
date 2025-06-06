# Databricks notebook source
# MAGIC %md
# MAGIC #### Last updated: 06/13/2024

# COMMAND ----------

# MAGIC %md
# MAGIC # FG: AD risk factors

# COMMAND ----------

# Normal analysis starts here
library(tidyverse)
library(lubridate)
raw_data_path = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/mingzhoufu@mednet.ucla.edu/AD_trajectory/data/"
load(file = paste0(raw_data_path, "FG_test/AD_risk/FG_AD_sample_demo_0610.rda"))
load(file = paste0(raw_data_path, "FG_test/AD_risk/case_visit_wide_0610.rda"))
print(dim(case_visit_wide)) # dim = (24473,1532)
load(file = paste0(raw_data_path, "FG_test/AD_risk/control_visit_wide_0610.rda"))
print(dim(control_visit_wide)) # dim = (24473,2565)

# COMMAND ----------

# MAGIC %md
# MAGIC ## Data prepare

# COMMAND ----------

# MAGIC %md
# MAGIC ### ICDs to test

# COMMAND ----------

# ICD in AD cases but not in controls
diff_icd = setdiff(names(case_visit_wide), names(control_visit_wide))
print(length(diff_icd)) # N = 41
# check if the frequency > 10 - no ICD freq > 10 so just drop these ICDs
for (i in 2:length(diff_icd)) {
  icd = diff_icd[i]
  freq_yes = as.numeric(table(case_visit_wide[[icd]])[2])
  if (freq_yes > 10) {
    print(paste0(icd, ": ", freq_yes))
  }
}

# COMMAND ----------

overlap_icd = intersect(names(case_visit_wide), names(control_visit_wide))
print(length(overlap_icd)) # N = 1491
control_visit_wide_short = control_visit_wide %>% select(all_of(overlap_icd))
case_visit_wide_short = case_visit_wide %>% select(all_of(overlap_icd))
combined_visit_wide = rbind(case_visit_wide_short, control_visit_wide_short)
print(dim(combined_visit_wide)) # dim = (48946,1491)

# COMMAND ----------

# MAGIC %md
# MAGIC ### Check for prevalence

# COMMAND ----------

# Calculate the sum of each column
column_sums = colSums(combined_visit_wide[,2:1491])
# Print the sum of each column (only keep icds with prevalence >=1%)
icd_prev1 = column_sums[column_sums >= 0.01*nrow(combined_visit_wide)] %>% names() 
print(length(icd_prev1)) # 586 
# finalize dataset
combined_visit_wide_short = combined_visit_wide %>% select("PatientID", all_of(icd_prev1))
print(dim(combined_visit_wide_short)) # dim = (48946,587)
save(combined_visit_wide_short, file = paste0(raw_data_path, "FG_test/AD_risk/combined_visit_wide_short_0610.rda"))

# COMMAND ----------

# MAGIC %md
# MAGIC ## FG modeling

# COMMAND ----------

# MAGIC %md
# MAGIC We ran this in sepearte files to save time. See `AD_FG_test_sep` for details.

# COMMAND ----------

# MAGIC %md
# MAGIC ### Clean results

# COMMAND ----------

library(tidyverse)
library(lubridate)
raw_data_path = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/mingzhoufu@mednet.ucla.edu/AD_trajectory/data/"
folder_path = paste0(raw_data_path, "FG_test/AD_risk/FG_result/tmp-results/") 
# List all .rda files from the folder
files = list.files(path = folder_path, pattern = "\\.rda$")
print(length(files)) # 8

# COMMAND ----------

for (i in 1:length(files)){
  load(paste0(folder_path, files[i]))
  if (i == 1) {
    FG_result_full = FG_result
  } else {
    FG_result_full = rbind(FG_result_full, FG_result)
  }
}
print(dim(FG_result_full)) # dim = (385,7)
save(FG_result_full, file = paste0(raw_data_path, "FG_test/AD_risk/FG_result/FG_result_AD.rda"))

# COMMAND ----------

load(file = paste0(raw_data_path, "FG_test/AD_risk/FG_result/FG_result_AD_0610.rda"))
head(FG_result_full)

# COMMAND ----------

AD_risk_tbl = FG_result_full %>% 
  mutate(coef = as.numeric(coef), `exp(coef)` = as.numeric(`exp(coef)`),
    z = as.numeric(z), `p-value` = as.numeric(`p-value`)) %>% filter(!is.na(coef)) %>%
  mutate(fdr_adj_p = p.adjust(`p-value`, method = "fdr")) %>% 
  filter(coef > 0 & fdr_adj_p <= 0.1) %>% arrange(desc(coef))
print(dim(AD_risk_tbl))
save(AD_risk_tbl, file = paste0(raw_data_path, "FG_test/AD_risk/FG_result/AD_risk_tbl_0610.rda"))
AD_risk_ICD = AD_risk_tbl %>% pull(exposure_ICD) %>% sort()
print(AD_risk_ICD)
head(AD_risk_tbl)

# COMMAND ----------

install.packages("icd.data")
library(icd.data)

# COMMAND ----------

# add description
icd10cm2016_3digit = icd10cm2016 %>% filter(nchar(code) == 3) %>% 
  mutate(code = as.character(code))
risk_AD_sig_desc = AD_risk_tbl %>% 
  left_join(icd10cm2016_3digit, by = c("exposure_ICD" = "three_digit")) %>% 
  mutate(HR = round(`exp(coef)`, 2)) %>%
  mutate(adj_p = case_when(
    fdf_adj_p < 0.001 ~ "<0.001",
    TRUE ~ as.character(round(fdf_adj_p, 3))
  )) %>%
  select(exposure_ICD, HR, adj_p, short_desc, chapter) 
risk_AD_sig_desc

# COMMAND ----------

# MAGIC %md
# MAGIC # FG: Other step pairs

# COMMAND ----------

load(file = paste0(raw_data_path, "patient_data/mod/diag_beforeAD_3digit_0610.rda"))
print(dim(diag_beforeAD_3digit)) # dim = (729347,4)
load(file = paste0(raw_data_path, "FG_test/AD_risk/case_visit_wide_0610.rda"))
print(dim(case_visit_wide)) # dim = (24473,1532)

# COMMAND ----------

# MAGIC %md
# MAGIC ## Finalize ICD list for testing

# COMMAND ----------

unique_icd_lst = unique(diag_beforeAD_3digit$ICD_3digit) 
print(length(unique_icd_lst)) # N = 1531
# need to overlap with the embedding distance matrix
embed_overlap = intersect(unique_icd_lst, embed_icd) %>% sort()
# also remove certain chapter codes
rm_codes = grep("^[OPQRVXY][0-9]+", embed_overlap, value = TRUE)
chap_rm_overlap = setdiff(embed_overlap, rm_codes) 
# -- deal with 1% prevalence --
# Calculate the sum of each column
column_sums = colSums(case_visit_wide[,2:1532])
# Print the sum of each column (only keep icds with prevalence >=1%)
icd_prev1 = column_sums[column_sums >= 0.01*nrow(case_visit_wide)] %>% names() 
print(length(icd_prev1)) # 459
icd_lst_final = intersect(chap_rm_overlap, icd_prev1)
print(length(icd_lst_final)) # N = 302
save(icd_lst_final, file = paste0(raw_data_path, "FG_test/other_risk/icd_lst_final_0610.rda"))

# COMMAND ----------

# clean patient records to restrict to the final icd list
diag_beforeAD_filter = diag_beforeAD_3digit %>% filter(ICD_3digit %in% icd_lst_final) %>% arrange(PatientID, EncounterDate)
print(dim(diag_beforeAD_filter)) # dim = (393601,4)
save(diag_beforeAD_filter, file = paste0(raw_data_path, "FG_test/other_risk/diag_beforeAD_filter_0610.rda"))

# COMMAND ----------

# MAGIC %md
# MAGIC ## Run FG models sepearately

# COMMAND ----------

# MAGIC %md
# MAGIC We ran this in sepearte files to save time. See `other_FG_test_sep` for details.

# COMMAND ----------

# MAGIC %md
# MAGIC ### Clean results

# COMMAND ----------

outcome_icd = "S90"
folder_path = paste0(raw_data_path, "FG_test/other_risk/FG_result/", outcome_icd, "/") 
# List all .rda files from the folder
files = list.files(path = folder_path, pattern = "\\.rda$")
print(length(files)) # 300
# files
for (i in 1:length(files)){
  load(paste0(folder_path, files[i]))
  if (i == 1) {
    FG_result = coefficient
  } else {
    FG_result = rbind(FG_result, coefficient)
  }
}
print(dim(FG_result))
print(head(FG_result))
save(FG_result, file = paste0(raw_data_path, "FG_test/other_risk/FG_result/", outcome_icd, ".rda"))

# COMMAND ----------

load(file = paste0(raw_data_path, "FG_test/other_risk/FG_result/FG_result_update.rda"))
other_risk_tbl = FG_result_update %>% 
  mutate(coef = as.numeric(coef), `exp(coef)` = as.numeric(`exp(coef)`),
    z = as.numeric(z), `p-value` = as.numeric(`p-value`)) %>% filter(!is.na(coef)) %>%
  mutate(fdr_adj_p = p.adjust(`p-value`, method = "fdr")) %>% 
  filter(coef > 0 & fdr_adj_p <= 0.1) %>% arrange(desc(coef))
print(dim(other_risk_tbl)) # dim = (32654,8)
save(other_risk_tbl, file = paste0(raw_data_path, "FG_test/other_risk/FG_result/other_risk_tbl_0610.rda"))

# COMMAND ----------

# MAGIC %md
# MAGIC ## Combine all FG results

# COMMAND ----------

load(file = paste0(raw_data_path, "FG_test/AD_risk/FG_result/AD_risk_tbl_0610.rda"))
# combine all
risk_pair_sig = rbind(AD_risk_tbl, other_risk_tbl)
print(dim(risk_pair_sig)) # dim = (32670,8)
save(risk_pair_sig, file = paste0(raw_data_path, "FG_test/risk_pair_sig_0610.rda"))

# COMMAND ----------

head(risk_pair_sig)

# COMMAND ----------

load(file = paste0(raw_data_path, "FG_test/risk_pair_sig_0610.rda"))

# COMMAND ----------

risk_pair_sig %>% filter(exposure_ICD == "F41", outcome_ICD == "F32")

# COMMAND ----------

risk_pair_sig %>% filter(exposure_ICD == "F32", outcome_ICD == "F41")

# COMMAND ----------

