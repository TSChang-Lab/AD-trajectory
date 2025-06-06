# Databricks notebook source
# MAGIC %md
# MAGIC #### Last updated: 03/28/2025

# COMMAND ----------

# MAGIC %md
# MAGIC # All UC data warehouse patients

# COMMAND ----------

# MAGIC %md
# MAGIC ## Clean demographic information

# COMMAND ----------

library(SparkR)
# All UC data warehouse patient information (run once)
all_pat_info = read.df("s3://uchdw-501868445017-us-west-2-prod-databricks-user-files/ucla_mingzhoufu/all_pat_info_0328", "csv", header = TRUE)
all_pat_info_df = collect(all_pat_info)
save(all_pat_info_df, file = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/mingzhoufu@mednet.ucla.edu/AD_trajectory/Review_Added/data/patient_data/raw/all_pat_info_df_0328.rda")

# COMMAND ----------

# Normal analysis starts here
library(tidyverse)
library(lubridate)
raw_data_path = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/mingzhoufu@mednet.ucla.edu/AD_trajectory/Review_Added/data/"
load(file = paste0(raw_data_path, "patient_data/raw/all_pat_info_df_0328.rda"))
print(dim(all_pat_info_df)) # dim = (9793287,7)
print(length(unique(all_pat_info_df$person_id))) # N = 9793287
head(all_pat_info_df)

# COMMAND ----------

# Clean patient demographics
pat_full_demo = all_pat_info_df %>% 
  mutate(BirthDate = as.Date(paste0(year_of_birth, '-06-30')), DeathDate = as.Date(death_date)) %>% 
  mutate(female = case_when(
    gender == "FEMALE" ~ 1,
    gender == "MALE" ~ 0
  )) %>%
  mutate(race_ethnicity = case_when(
    ethnicity == "Hispanic or Latino" ~ "Hispanic",
    race == "White" & ethnicity == "Not Hispanic or Latino" ~ "NH-White",
    race == "Black or African American" ~ "Black",
    race == "Asian" ~ "Asian",
    race == "American Indian or Alaska Native" ~ "Amerian-Indian",
    race %in% c("Native Hawaiian or Other Pacific Islander", "Multirace", "Other Race") ~ "Others",
    TRUE ~ "Unknown"
  )) %>% rename("PatientID" = "person_id") %>%
  mutate(PatientLivingStatus = if_else(!is.na(DeathDate), "Deceased", "Alive")) %>%
  select(PatientID, female, race_ethnicity, PatientLivingStatus, location_source_value, BirthDate, DeathDate) %>%
  filter(!is.na(female) & race_ethnicity != "Unknown")
print(dim(pat_full_demo)) # dim = (6741476,7)
save(pat_full_demo, file = paste0(raw_data_path, "patient_data/mod/pat_full_demo_0328.rda"))

# COMMAND ----------

head(pat_full_demo)

# COMMAND ----------

# output to file
pat_clean_id = pat_full_demo %>% select(PatientID) %>% unique()
print(dim(pat_clean_id)) # dim = (6741476,1)
# output to csv
write.table(pat_clean_id, paste0(raw_data_path, "patient_data/mod/pat_clean_id_0328.csv"), sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE)

# COMMAND ----------

# MAGIC %md
# MAGIC ## Clean EHR features

# COMMAND ----------

library(SparkR)
# All UC data warehouse patient information (run once)
all_pat_visit_info = read.df("s3://uchdw-501868445017-us-west-2-prod-databricks-user-files/ucla_mingzhoufu/all_pat_visit_0328", "csv", header = TRUE)
all_pat_visit_info_df = SparkR::collect(all_pat_visit_info)
save(all_pat_visit_info_df, file = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/mingzhoufu@mednet.ucla.edu/AD_trajectory/Review_Added/data/patient_data/raw/all_pat_visit_info_df_0328.rda")

# COMMAND ----------

# Normal analysis starts here
library(tidyverse)
library(lubridate)
raw_data_path = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/mingzhoufu@mednet.ucla.edu/AD_trajectory/Review_Added/data/"
load(file = paste0(raw_data_path, "patient_data/raw/all_pat_visit_info_df_0328.rda"))
print(dim(all_pat_visit_info_df)) # dim = (5667122,5)
load(file = paste0(raw_data_path, "patient_data/mod/pat_full_demo_0328.rda"))
print(dim(pat_full_demo)) # dim = (6741476,7)

# COMMAND ----------

# MAGIC %md
# MAGIC ## Select eligible UC patients

# COMMAND ----------

uc_pat_full = all_pat_visit_info_df %>% rename("PatientID" = "person_id") %>%
  mutate(date_first_visit = as.Date(date_first_visit), date_last_visit = as.Date(date_last_visit)) %>% 
  inner_join(pat_full_demo) %>% mutate(n_encounter = as.numeric(n_encounter), n_icd_3digit = as.numeric(n_icd_3digit)) %>%
  mutate(age_last_visit = as.numeric(date_last_visit - BirthDate)/365.25,
        record_length = as.numeric(date_last_visit - date_first_visit)/365.25, 
        enc_per_yr = n_encounter/record_length) %>% 
  select(PatientID, female, race_ethnicity, PatientLivingStatus, location_source_value, BirthDate, DeathDate,
        date_first_visit, date_last_visit, n_encounter, n_icd_3digit, age_last_visit, record_length, enc_per_yr)
print(dim(uc_pat_full)) # dim = (5667122,14)

# COMMAND ----------

uc_pat_eligible = uc_pat_full %>% filter(age_last_visit >= 65 & age_last_visit < 90) %>%
  filter(record_length <= age_last_visit & n_encounter >= 2 & enc_per_yr >= 1)
print(dim(uc_pat_eligible)) # dim = (970348,14)
save(uc_pat_eligible, file = paste0(raw_data_path, "patient_data/final/uc_pat_eligible_0328.rda"))

# COMMAND ----------

head(uc_pat_eligible)

# COMMAND ----------

# output IDs to file
pat_eligible_id = uc_pat_eligible %>% select(PatientID) %>% unique()
print(dim(pat_eligible_id)) # dim = (970348,1)
# output to csv
write.table(pat_eligible_id, paste0(raw_data_path, "patient_data/mod/pat_eligible_id_0328.csv"), sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE)

# COMMAND ----------

# MAGIC %md
# MAGIC # AD patient visit cleaning

# COMMAND ----------

library(SparkR)
# All AD patient information (run once)
AD_pat_visit = read.df("s3://uchdw-501868445017-us-west-2-prod-databricks-user-files/ucla_mingzhoufu/AD_encounters_0328", "csv", header = TRUE)
AD_pat_visit_df = SparkR::collect(AD_pat_visit)
save(AD_pat_visit_df, file = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/mingzhoufu@mednet.ucla.edu/AD_trajectory/Review_Added/data/patient_data/raw/AD_pat_visit_df_0328.rda")

# COMMAND ----------

# Normal analysis starts here
library(tidyverse)
library(lubridate)
'%!in%' = function(x,y)!{'%in%'(x,y)}
raw_data_path = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/mingzhoufu@mednet.ucla.edu/AD_trajectory/Review_Added/data/"
load(file = paste0(raw_data_path, "patient_data/raw/AD_pat_visit_df_0328.rda"))
print(dim(AD_pat_visit_df)) # dim = (2021588,3)
print(length(unique(AD_pat_visit_df$person_id))) # 31914 | 33899
head(AD_pat_visit_df)

# COMMAND ----------

load(file = paste0(raw_data_path, "patient_data/final/uc_pat_eligible_0328.rda"))
# clean the data
AD_encounters_raw = AD_pat_visit_df %>%
  mutate(EncounterDate = as.Date(condition_start_date)) %>% rename("PatientID" = "person_id", "ICD" = "icd") %>%
  select(PatientID, EncounterDate, ICD) %>% unique()
print(dim(AD_encounters_raw)) # dim = (2021588,3)
# define AD date
AD_icd = c("G30", "G30.0", "G30.1", "G30.8", "G30.9")
AD_demographic = AD_encounters_raw %>% filter(ICD %in% AD_icd) %>% group_by(PatientID) %>%
  summarise(AD_Date = as.Date(min(EncounterDate))) %>% ungroup() %>% inner_join(uc_pat_eligible) %>% unique()
print(dim(AD_demographic)) # dim = (25977,15)

# COMMAND ----------

# there are 3 patients with different death date, we remove them for consistency
rm_pat = AD_demographic %>% group_by(PatientID) %>% mutate(n = n()) %>% filter(n > 1) %>% pull(PatientID) %>% unique()
print(length(rm_pat)) # N = 0
AD_demographic = AD_demographic %>% filter(PatientID %!in% rm_pat) %>%
  filter(AD_Date <= as.Date("2024-05-10"))
print(dim(AD_demographic)) # dim = (24824,15)
save(AD_demographic, file = paste0(raw_data_path, "patient_data/mod/AD_demographic_0328.rda"))

# COMMAND ----------

head(AD_demographic)

# COMMAND ----------

AD_encounters_raw_add = AD_encounters_raw %>% left_join(AD_demographic) %>% 
  select(PatientID, EncounterDate, ICD, AD_Date) %>% filter(EncounterDate == AD_Date) %>%
  select(PatientID, ICD) %>% unique() %>% mutate(status = 1) %>%
  pivot_wider(names_from = ICD, values_from = status, values_fill = 0)
print(dim(AD_encounters_raw_add)) # (24824,7216)
head(AD_encounters_raw_add)

# COMMAND ----------

# Calculate the sum of each column
column_sums = colSums(AD_encounters_raw_add[,2:7216])/length(unique(AD_encounters_raw_add$PatientID))
# filter out co-occur with AD (>0.5)
column_sums[column_sums > 0.5]

# COMMAND ----------

F028_encounters = AD_encounters_raw %>% left_join(AD_demographic) %>% 
  select(PatientID, EncounterDate, ICD, AD_Date) %>% filter(substr(ICD, 1, 5) == "F02.8") %>% 
  filter(EncounterDate <= AD_Date) %>%
  mutate(time_diff = as.numeric(AD_Date - EncounterDate)) %>% filter(!is.na(AD_Date))
print(dim(F028_encounters)) # dim = (24148,5)
print(length(unique(F028_encounters$PatientID)))  # 22686
head(F028_encounters)

# COMMAND ----------

summary(F028_encounters$time_diff)

# COMMAND ----------

sum(F028_encounters$time_diff >= 90) # 1204

# COMMAND ----------

# MAGIC %md
# MAGIC Given F02.8 usually (72% of all time) happens at the same time with AD (G30), and only 5% (1204/22686) of time happened three months or more before AD, we decided to exclude this diagnosis (within 3 months with AD) in the following analysis.

# COMMAND ----------

# MAGIC %md
# MAGIC ## AD cases vs. control sample prepare

# COMMAND ----------

# MAGIC %md
# MAGIC ### AD cases and non-AD controls cleaning

# COMMAND ----------

# AD cases
# Filter conditions no later than AD and recalculate EHR features
AD_date = AD_demographic %>% select(PatientID, AD_Date, date_first_visit) %>% unique()
AD_encounters_filter = AD_encounters_raw %>% inner_join(AD_date) %>% filter(EncounterDate <= AD_Date) %>% 
  mutate(rmv = case_when(
    substr(ICD, 1, 5) == "F02.8" & (AD_Date - EncounterDate <= 90) ~ "yes",
    TRUE ~ "no"
  )) %>% filter(rmv == "no") %>% select(-rmv)
print(dim(AD_encounters_filter)) # dim = (996418,5)
save(AD_encounters_filter, file = paste0(raw_data_path, "patient_data/mod/AD_encounters_filter_0328.rda"))
# calculate new EHR features
AD_pat_EHR_feature = AD_encounters_filter %>% group_by(PatientID) %>%
  summarise(n_encounter = length(unique(EncounterDate)), follow_up_time = as.numeric(AD_Date - date_first_visit)/365.25, 
        enc_per_yr = n_encounter/follow_up_time) %>% ungroup() %>% unique()
print(dim(AD_pat_EHR_feature)) # dim = (24824,4)
# finalize AD cases info
case_demo = AD_demographic %>% 
  select(PatientID, female, race_ethnicity, PatientLivingStatus, location_source_value, age_last_visit) %>% inner_join(AD_pat_EHR_feature)
print(dim(case_demo)) # dim = (24824,9)

# COMMAND ----------

# non-AD controls (further restrict to longer records and higher density)
control_demo = uc_pat_eligible %>% filter(PatientID %!in% unique(AD_pat_visit_df$person_id)) %>%
  filter(record_length >= 5 & n_encounter >= 2 & enc_per_yr >= 2) %>% unique() %>%
  rename(follow_up_time = record_length) %>% select(all_of(names(case_demo)))
print(dim(control_demo)) # dim = (313971,9)

# COMMAND ----------

# MAGIC %md
# MAGIC ### Matching

# COMMAND ----------

install.packages("MatchIt")
library(MatchIt)

# COMMAND ----------

case_demo$case = 1
control_demo$case = 0
full_demo = rbind(case_demo, control_demo)
# perform matching
match_obj = matchit(case ~ female + race_ethnicity + location_source_value + age_last_visit,
  data = full_demo, method = "nearest", distance = "glm", ratio = 3, replace = FALSE)
summary(match_obj)

# COMMAND ----------

# finalize the data
matched_AD_control = match.data(match_obj) %>% mutate(group = if_else(case == 1, "case", "control")) %>%
  select(all_of(names(case_demo)), group) %>% select(-case)
print(dim(matched_AD_control)) # dim = (99296,10)
head(matched_AD_control)
save(matched_AD_control, file = paste0(raw_data_path, "patient_data/final/matched_AD_control_0328.rda"))

# COMMAND ----------

install.packages("arsenal")
library(arsenal)

# COMMAND ----------

table_one = tableby(group ~ age_last_visit + as.factor(female) + race_ethnicity + location_source_value, 
  data = matched_AD_control, numeric.stats = "meansd") 
summary(table_one, title = "Matched Data")

# COMMAND ----------

# output to file
matched_control_id = matched_AD_control %>% filter(group == "control") %>% select(PatientID) %>% unique()
print(dim(matched_control_id)) # dim = (74472,1)
# output to csv
write.table(matched_control_id, paste0(raw_data_path, "patient_data/mod/matched_control_id_0328.csv"), sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE)

# COMMAND ----------

library(SparkR)
# All matched control patient information (run once)
control_pat_visit = read.df("s3://uchdw-501868445017-us-west-2-prod-databricks-user-files/ucla_mingzhoufu/control_encounters_0328", "csv", header = TRUE)
control_pat_visit_df = SparkR::collect(control_pat_visit)
save(control_pat_visit_df, file = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/mingzhoufu@mednet.ucla.edu/AD_trajectory/Review_Added/data/patient_data/raw/control_pat_visit_df_0328.rda")

# COMMAND ----------

load(file = paste0(raw_data_path, "patient_data/raw/control_pat_visit_df_0328.rda"))
# check F02.8 frequency in the control
control_F028 = control_pat_visit_df %>% filter(substr(icd, 1, 5) == "F02.8")
print(dim(control_F028)) # dim = (1717,3)

# COMMAND ----------

# MAGIC %md
# MAGIC ### FG sample prepare (AD & non-AD)

# COMMAND ----------

# Normal analysis starts here
library(tidyverse)
library(lubridate)
'%!in%' = function(x,y)!{'%in%'(x,y)}
raw_data_path = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/mingzhoufu@mednet.ucla.edu/AD_trajectory/Review_Added/data/"
load(file = paste0(raw_data_path, "patient_data/final/matched_AD_control_0328.rda"))
print(dim(matched_AD_control)) # dim = (99296,10)
load(file = paste0(raw_data_path, "patient_data/final/uc_pat_eligible_0328.rda"))

# COMMAND ----------

age_first_visit_add = uc_pat_eligible %>% 
  mutate(age_first_visit = as.numeric(date_first_visit - BirthDate)/365.25) %>%
  select(PatientID, age_first_visit) %>% unique()
print(dim(age_first_visit_add)) # dim = (970348, 2)
save(age_first_visit_add, file = paste0(raw_data_path, "patient_data/mod/age_first_visit_add_0328.rda"))
# get full covariates
matched_AD_control_demo = matched_AD_control %>% left_join(age_first_visit_add) %>%
  mutate(outcome_status = case_when(
    group == "control" & PatientLivingStatus == "Deceased" ~ 2,
    group == "case" ~ 1,
    TRUE ~ 0
  )) %>% select(PatientID, age_first_visit, female, race_ethnicity, follow_up_time, outcome_status) %>% unique()
print(dim(matched_AD_control_demo)) # dim = (99296, 6)
print(table(matched_AD_control_demo$outcome_status))
save(matched_AD_control_demo, file = paste0(raw_data_path, "FG_test/AD_risk/matched_AD_control_demo_0328.rda"))

# COMMAND ----------

# FG modeling need follow_up_time > 0
FG_AD_sample_demo = matched_AD_control_demo %>% filter(follow_up_time > 0)
print(table(FG_AD_sample_demo$outcome_status))
save(FG_AD_sample_demo, file = paste0(raw_data_path, "FG_test/AD_risk/FG_AD_sample_demo_0328.rda"))

# COMMAND ----------

# MAGIC %md
# MAGIC ## Exposure status prepare

# COMMAND ----------

# MAGIC %md
# MAGIC ### AD patient visits

# COMMAND ----------

load(file = paste0(raw_data_path, "patient_data/mod/AD_encounters_filter_0328.rda"))
diag_beforeAD_3digit = AD_encounters_filter %>% mutate(ICD_3digit = substr(ICD, 1, 3)) %>%
  select(-c(ICD, date_first_visit)) %>% unique() %>% group_by(PatientID, ICD_3digit) %>%
  arrange(EncounterDate) %>% slice(1) %>% ungroup() %>% filter(EncounterDate <= AD_Date)
print(dim(diag_beforeAD_3digit)) # dim = (743946,4)
save(diag_beforeAD_3digit, file = paste0(raw_data_path, "patient_data/mod/diag_beforeAD_3digit_0328.rda"))
# combine all together to make a big dummy dataset
case_visit_wide = diag_beforeAD_3digit %>%
  select(PatientID, ICD_3digit) %>% mutate(status = 1) %>% unique() %>%
  pivot_wider(names_from = ICD_3digit, values_from = status, values_fill = 0) %>% 
  mutate_all(~ifelse(is.na(.), 0, .)) 
print(dim(case_visit_wide)) # dim = (24824,1547)
save(case_visit_wide, file = paste0(raw_data_path, "FG_test/AD_risk/case_visit_wide_0328.rda"))

# COMMAND ----------

# MAGIC %md
# MAGIC ### Control patient visits

# COMMAND ----------

load(file = paste0(raw_data_path, "patient_data/raw/control_pat_visit_df_0328.rda"))
print(dim(control_pat_visit_df)) # dim = (9585405,3)
# clean the data
control_encounters_raw = control_pat_visit_df %>%
  mutate(EncounterDate = as.Date(condition_start_date)) %>% rename("PatientID" = "person_id", "ICD" = "icd") %>%
  select(PatientID, EncounterDate, ICD) %>% unique()
print(dim(control_encounters_raw)) # dim = (9585405,3)
# cut to 3-digit level ICD
diag_control_3digit = control_encounters_raw %>% mutate(ICD_3digit = substr(ICD, 1, 3)) %>%
  select(-c(ICD)) %>% unique() %>% group_by(PatientID, ICD_3digit) %>%
  arrange(EncounterDate) %>% slice(1) %>% ungroup() 
print(dim(diag_control_3digit)) # dim = (6673331,3)
save(diag_control_3digit, file = paste0(raw_data_path, "patient_data/mod/diag_control_3digit_0328.rda"))

# COMMAND ----------

# combine all together to make a big dummy dataset
control_visit_wide = diag_control_3digit %>%
  select(PatientID, ICD_3digit) %>% mutate(status = 1) %>% unique() %>%
  pivot_wider(names_from = ICD_3digit, values_from = status, values_fill = 0) %>% 
  mutate_all(~ifelse(is.na(.), 0, .)) 
print(dim(control_visit_wide)) # dim = (74472,2763)
save(control_visit_wide, file = paste0(raw_data_path, "FG_test/AD_risk/control_visit_wide_0328.rda"))

# COMMAND ----------

head(matched_AD_control)

# COMMAND ----------

# MAGIC %md
# MAGIC # Sample comparison

# COMMAND ----------

library(tidyverse)
'%!in%' = function(x,y)!{'%in%'(x,y)}
raw_data_path = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/mingzhoufu@mednet.ucla.edu/AD_trajectory/data/"
load(file = paste0(raw_data_path, "patient_data/mod/AD_demographic_0610.rda"))
print(dim(AD_demographic)) # dim = (24473,15)
# Load in clustered patients
load(file = paste0(raw_data_path, "patient_data/final/pat_traj_final_cluster_0622.rda"))
print(dim(pat_traj_final_cluster)) # dim = (13409,14)

# COMMAND ----------

pat_info_UC = pat_traj_final_cluster %>% left_join(AD_demographic) %>%
  mutate(age_AD = as.numeric(AD_Date - BirthDate)/365.25) %>%
  mutate(AD_to_last = as.numeric(date_last_visit - AD_Date)/365.25) %>%
  mutate(AD_to_death = as.numeric(DeathDate - AD_Date)/365.25) %>%
  mutate(first_visit_to_AD = as.numeric(AD_Date - date_first_visit)/365.25) %>%
  select(female, race_ethnicity, PatientLivingStatus, n_icd_3digit, n_encounter, record_length, enc_per_yr, 
         age_last_visit, first_visit_to_AD, age_AD, AD_to_last, AD_to_death) %>% unique() %>%
  mutate(female = as.factor(female))
print(dim(pat_info_UC)) # dim = (5762,12)

# COMMAND ----------

cont_summary = pat_info_UC %>% 
  select(where(is.numeric)) %>% 
  summarise(across(
    .fns = list(
      median = ~ median(.x, na.rm = TRUE),
      Q1     = ~ quantile(.x, .25, na.rm = TRUE),
      Q3     = ~ quantile(.x, .75, na.rm = TRUE)
    )
  ), .groups = "drop") %>%
  # pivot longer for nicer display
  pivot_longer(
    cols      = everything(), 
    names_to  = c("variable","stat"),
    # regex: capture everything up to the last underscore, then the suffix
    names_pattern = "(.*)_(median|Q1|Q3)$"
  ) %>%
  pivot_wider(names_from = stat, values_from = value)
cont_summary

# COMMAND ----------

cat_summary = pat_info_UC %>%
  select(where(~ is.character(.x) || is.factor(.x))) %>%
  summarise(across(everything(), ~ list(table(.x)))) %>%
  tidyr::pivot_longer(everything(), 
                      names_to = "variable", 
                      values_to = "tbl") %>%
  rowwise() %>%
  mutate(
    counts = list(as.integer(tbl)),
    levels = list(names(tbl)),
    perc    = list( round(100 * counts / sum(counts), 1) )
  ) %>%
  select(variable, levels, counts, perc)
cat_summary_flat = cat_summary %>%
  unnest(cols = c(levels, counts, perc)) %>%
  rename(
    level   = levels,
    count   = counts,
    percent = perc
  )

print(cat_summary_flat)

# COMMAND ----------

