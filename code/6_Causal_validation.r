# Databricks notebook source
# MAGIC %md
# MAGIC #### Last updated: 07/02/2024

# COMMAND ----------

library(tidyverse)
'%!in%' = function(x,y)!{'%in%'(x,y)}
raw_data_path = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/mingzhoufu@mednet.ucla.edu/AD_trajectory/data/"

# COMMAND ----------

# MAGIC %md
# MAGIC # Causal structural learning

# COMMAND ----------

# MAGIC %md
# MAGIC ## AD links (AD + controls)

# COMMAND ----------

# MAGIC %md
# MAGIC ### Data preparation

# COMMAND ----------

load(file = paste0(raw_data_path, "FG_test/AD_risk/combined_visit_wide_short_0610.rda"))
case_control_prep = combined_visit_wide_short %>% select(-PatientID)
print(dim(case_control_prep)) # dim = (48946,586)

# COMMAND ----------

# Calculate the sum of each column
column_sums = colSums(case_control_prep)
# Print the sum of each column (only keep icds with prevalence >=1%)
select_icds = column_sums[column_sums >= 0.01*nrow(case_control_prep)] %>% names() 
print(length(select_icds)) # 586
load(file = paste0(raw_data_path, "upload_ucla/embed_sim_cal.rda"))
embed_icd = intersect(select_icds, colnames(embed_sim_cal)) %>% sort()
print(length(embed_icd)) # 462
# also remove certain chapter codes
rm_codes = grep("^[OPQRVXY][0-9]+", embed_icd, value = TRUE)
icd_lst_filter = setdiff(embed_icd, rm_codes) 
print(length(icd_lst_filter)) # N = 388
causal_pat_short = case_control_prep %>% select(all_of(icd_lst_filter))
print(dim(causal_pat_short)) # dim = (48946,388)
causal_pat_short$G30 = c(rep(1, 24473), rep(0, 24473))

# COMMAND ----------

# MAGIC %md
# MAGIC ### Split into 10 folds

# COMMAND ----------

case_pat_short = causal_pat_short[1:24473,]
# Calculate the number of rows in each subset
subset_size = ceiling(nrow(case_pat_short) / 10)
# Create an index vector to assign each row to a subset
index = rep(1:10, each = subset_size)
# Adjust the index vector to handle cases where the dataset size is not divisible by 10
index = index[1:nrow(case_pat_short)]
# Split the dataset into ten subsets based on the index vector
subset_list_case = split(case_pat_short, index)

# COMMAND ----------

control_pat_short = causal_pat_short[24474:48946,]
# Calculate the number of rows in each subset
subset_size = ceiling(nrow(control_pat_short) / 10)
# Create an index vector to assign each row to a subset
index = rep(1:10, each = subset_size)
# Adjust the index vector to handle cases where the dataset size is not divisible by 10
index = index[1:nrow(control_pat_short)]
# Split the dataset into ten subsets based on the index vector
subset_list_control = split(control_pat_short, index)

# COMMAND ----------

# Initialize an empty list to store combined data frames
combined_datasets = list()
# Iterate through each data frame in the list
for (i in 1:10) {
  # Exclude the i-th data frame
  combined_cases = do.call(rbind, subset_list_case[-i])
  combined_controls = do.call(rbind, subset_list_control[-i])
  # Store the combined data frame in the list
  combined_datasets[[i]] = rbind(combined_cases, combined_controls)
}

# COMMAND ----------

save(combined_datasets, file = paste0(raw_data_path, "causal/AD_links/combined_datasets_CV_0703.rda"))

# COMMAND ----------

# MAGIC %md
# MAGIC ### Run GES

# COMMAND ----------

library(remotes)
install_version("BiocManager", version = "1.30.22")
BiocManager::install(version = "3.18")
BiocManager::install("BiocGenerics")
BiocManager::install("graph")
BiocManager::install("RBGL")
install.packages("pcalg")
library(pcalg)

# COMMAND ----------

# MAGIC %md
# MAGIC We run this in the background with 10-fold cross validation. See details `AD_GES`.

# COMMAND ----------

# MAGIC %md
# MAGIC ### Clean causal results

# COMMAND ----------

for (n in 1:10) {
  load(file = paste0(raw_data_path, "causal/AD_links/ges.fit_", n, ".rda"))
  all_nodes = ges.fit$essgraph$`.->.nodes`
  arc_array = ges.fit$essgraph$`.->.in.edges`
  for (i in 1:length(arc_array)) {
    test_node = all_nodes[i]
    node_arc = arc_array[[test_node]]
    causal_learned = cbind(all_nodes[node_arc], rep(test_node, length(node_arc))) %>% 
      as.data.frame() 
    if (i == 1) {causal_learned_df = causal_learned}
    else {causal_learned_df = rbind(causal_learned_df, causal_learned)}
  }
  names(causal_learned_df) = c("start", "end")
  causal_learned_df$n_iter = n
  if (n == 1) {causal_full = causal_learned_df}
  else {causal_full = rbind(causal_full, causal_learned_df)}
}
print(dim(causal_full)) # dim = (87365,3)
save(causal_full, file = paste0(raw_data_path, "causal/AD_links/causal_full_AD_0709.rda"))

# COMMAND ----------

load(file = paste0(raw_data_path, "causal/AD_links/causal_full_AD_0709.rda"))

# COMMAND ----------

summary = causal_full %>% group_by(start, end) %>% summarize(n_time = n()) %>% ungroup()
table(summary$n_time)

# COMMAND ----------

# clean up the results (select links with >5 times)
causal_filtered_AD = causal_full %>% group_by(start, end) %>% 
  summarize(n_time = n()) %>% ungroup() %>% filter(n_time >= 1)
causal_filtered_AD = causal_filtered_AD %>% select(start, end) %>% unique() %>% mutate(causal = 1)
print(dim(causal_filtered_AD)) # dim = (18907,3)
save(causal_filtered_AD, file = paste0(raw_data_path, "causal/AD_links/causal_filtered_AD_0709.rda"))

# COMMAND ----------

causal_full %>% filter(end == "G30") %>% group_by(start) %>% summarise(n = n_distinct(n_iter)) %>% unique() %>% arrange(desc(n)) %>% filter(substr(start, 1, 1) %in% c("F", "G" ,"I")) 

# COMMAND ----------

load(file = paste0(raw_data_path, "clusters/filter_edges_1.rda"))
edge_df_depression = filter_edges %>% select(start, end) %>% 
  mutate(cluster = "Depression", asso_raw = 1)
load(file = paste0(raw_data_path, "clusters/filter_edges_2.rda"))
edge_df_vascular = filter_edges %>% select(start, end) %>% 
  mutate(cluster = "Vascular", asso_raw = 1)
load(file = paste0(raw_data_path, "clusters/filter_edges_3.rda"))
edge_df_encephalopathy = filter_edges %>% select(start, end) %>% 
  mutate(cluster = "Encephalopathy", asso_raw = 1)
load(file = paste0(raw_data_path, "clusters/filter_edges_4.rda"))
edge_df_mci = filter_edges %>% select(start, end) %>% 
  mutate(cluster = "MCI", asso_raw = 1)

# COMMAND ----------

edge_df_final = rbind(edge_df_depression, edge_df_vascular, edge_df_encephalopathy, edge_df_mci) %>% 
  mutate(value = 1) %>% unique() %>% pivot_wider(names_from = cluster,
    values_from = value,
    values_fill = list(value = 0)) %>% arrange(start, end)
print(dim(edge_df_final)) # dim = (130,7)
head(edge_df_final)

# COMMAND ----------

causal_merge = edge_df_final %>% left_join(causal_filtered_AD, by = c("start", "end")) %>% rename("causal_AD" = "causal")
print(dim(causal_merge)) # dim = (130,8)
head(causal_merge)

# COMMAND ----------

download_path = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/output/2024_09_16/"
save(causal_merge, file = paste0(download_path, "causal_merge_AD.rda"))

# COMMAND ----------

# MAGIC %md
# MAGIC ## Other links

# COMMAND ----------

# MAGIC %md
# MAGIC ### Data preparation

# COMMAND ----------

load(file = paste0(raw_data_path, "FG_test/AD_risk/case_visit_wide_0610.rda"))
case_prep = case_visit_wide %>% select(-PatientID)
print(dim(case_prep)) # dim = (24473,1531)

# COMMAND ----------

# Calculate the sum of each column
column_sums = colSums(case_prep)
# Print the sum of each column (only keep icds with prevalence >=1%)
select_icds = column_sums[column_sums >= 0.01*nrow(case_prep)] %>% names() 
print(length(select_icds)) # 459
causal_pat_short = case_prep %>% select(all_of(select_icds))
print(dim(causal_pat_short)) # dim = (24473,459)

# COMMAND ----------

# MAGIC %md
# MAGIC ### Split into 10 folds

# COMMAND ----------

# Calculate the number of rows in each subset
subset_size = ceiling(nrow(causal_pat_short) / 10)
# Create an index vector to assign each row to a subset
index = rep(1:10, each = subset_size)
# Adjust the index vector to handle cases where the dataset size is not divisible by 10
index = index[1:nrow(causal_pat_short)]
# Split the dataset into ten subsets based on the index vector
subset_list_case = split(causal_pat_short, index)
# Initialize an empty list to store combined data frames
combined_datasets = list()
# Iterate through each data frame in the list
for (i in 1:10) {
  # Exclude the i-th data frame
  combined_sample = do.call(rbind, subset_list_case[-i])
  # Store the combined data frame in the list
  combined_datasets[[i]] = combined_sample
}

# COMMAND ----------

save(combined_datasets, file = paste0(raw_data_path, "causal/other_links/combined_datasets_CV_ADonly_0627.rda"))

# COMMAND ----------

# MAGIC %md
# MAGIC ### Run GES

# COMMAND ----------

# MAGIC %md
# MAGIC We run this in the background with 10-fold cross validation. See details `other_GES`.

# COMMAND ----------

# MAGIC %md
# MAGIC ### Clean causal results

# COMMAND ----------

for (n in 1:10) {
  load(file = paste0(raw_data_path, "causal/other_links/ges.fit_", n, ".rda"))
  all_nodes = ges.fit$essgraph$`.->.nodes`
  arc_array = ges.fit$essgraph$`.->.in.edges`
  for (i in 1:length(arc_array)) {
    test_node = all_nodes[i]
    node_arc = arc_array[[test_node]]
    causal_learned = cbind(all_nodes[node_arc], rep(test_node, length(node_arc))) %>% 
      as.data.frame() 
    if (i == 1) {causal_learned_df = causal_learned}
    else {causal_learned_df = rbind(causal_learned_df, causal_learned)}
  }
  names(causal_learned_df) = c("start", "end")
  causal_learned_df$n_iter = n
  if (n == 1) {causal_full = causal_learned_df}
  else {causal_full = rbind(causal_full, causal_learned_df)}
}
print(dim(causal_full)) # dim = (83560,3)
save(causal_full, file = paste0(raw_data_path, "causal/other_links/causal_full_other_0702.rda"))

# COMMAND ----------

load(file = paste0(raw_data_path, "causal/other_links/causal_full_other_0702.rda"))

# COMMAND ----------

summary = causal_full %>% group_by(start, end) %>% summarize(n_time = n()) %>% ungroup()
table(summary$n_time)

# COMMAND ----------

# clean up the results (select links with >5 times)
causal_filtered_other = causal_full %>% group_by(start, end) %>% 
  summarize(n_time = n()) %>% ungroup() %>% filter(n_time >= 5)
causal_filtered_other = causal_filtered_other %>% select(start, end) %>% unique() %>% mutate(causal = 1)
print(dim(causal_filtered_other)) # dim = (5516,3)
save(causal_filtered_other, file = paste0(raw_data_path, "causal/other_links/causal_filtered_other_0702.rda"))

# COMMAND ----------

load(file = paste0(raw_data_path, "causal/other_links/causal_filtered_other_0702.rda"))

# COMMAND ----------

causal_merge_final = causal_merge %>% left_join(causal_filtered_other, by = c("start", "end")) %>% 
  rename("causal_other" = "causal") %>% mutate(casual = case_when(
    end == "G30" ~ causal_AD,
    end != "G30" ~ causal_other
  )) %>% select(-c(causal_AD, causal_other)) %>%
  mutate(casual_ges = case_when(
    is.na(casual) ~ 0,
    !is.na(casual) ~ casual
  )) %>% select(-casual)
print(dim(causal_merge_final)) # dim = (130,8)
head(causal_merge_final)

# COMMAND ----------

download_path = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/output/2024_09_16/"
save(causal_merge_final, file = paste0(download_path, "causal_merge_final.rda"))

# COMMAND ----------

write.csv(causal_merge_final, file = paste0(download_path, "causal_merge_final.csv"), row.names = FALSE)

# COMMAND ----------

table(causal_merge_final$casual_ges, useNA = "ifany")

# COMMAND ----------

# raw causal% (26.2%)
34/(34+96)

# COMMAND ----------

causal_merge_final %>% filter(casual_ges == 1 & end == "G30")

# COMMAND ----------

causal_merge_final %>% filter(casual_ges == 1 & Depression == 1) # 10

# COMMAND ----------

print(10/dim(edge_df_depression)[1])

# COMMAND ----------

causal_merge_final %>% filter(casual_ges == 1 & Encephalopathy == 1) # 18

# COMMAND ----------

print(18/dim(edge_df_encephalopathy)[1])

# COMMAND ----------

causal_merge_final %>% filter(casual_ges == 1 & MCI == 1) # 11

# COMMAND ----------

print(11/dim(edge_df_mci)[1])

# COMMAND ----------

/Workspace/Users/mingzhoufu@mednet.ucla.edu/output/2024_09_16

# COMMAND ----------

causal_merge_final %>% filter(casual_ges == 1 & Vascular == 1) # 14

# COMMAND ----------

print(14/dim(edge_df_vascular)[1])

# COMMAND ----------



# COMMAND ----------

# MAGIC %md
# MAGIC # Trajectory as risk factors

# COMMAND ----------

install.packages("cmprsk")
library(cmprsk)

# COMMAND ----------

load(file = paste0(raw_data_path, "patient_data/mod/diag_beforeAD_3digit_0610.rda"))
case_encounters = diag_beforeAD_3digit %>% select(-AD_Date)
print(dim(case_encounters)) # dim = (729347,3)
load(file = paste0(raw_data_path, "patient_data/mod/diag_control_3digit_rm_0703.rda"))
control_encounters = diag_control_3digit_rm
print(dim(control_encounters)) # dim = (1508752,3)
full_encounters = rbind(case_encounters, control_encounters) %>% as.data.frame()
print(dim(full_encounters)) # dim = (2238099,3)
print(length(unique(full_encounters$PatientID))) # 48945
head(full_encounters)

# COMMAND ----------

load(file = paste0(raw_data_path, "FG_test/AD_risk/FG_AD_sample_demo_rm_0703.rda"))
head(FG_AD_sample_demo_rm)

# COMMAND ----------

# MAGIC %md
# MAGIC ### Depression cluster

# COMMAND ----------

load(file = paste0(raw_data_path, "clusters/risk_test_depression.rda"))
head(risk_test_depression)

# COMMAND ----------

# MAGIC %md
# MAGIC ### Vascular cluster

# COMMAND ----------

load(file = paste0(raw_data_path, "clusters/risk_test_vascular.rda"))
head(risk_test_vascular)

# COMMAND ----------

# MAGIC %md
# MAGIC ### Encephalopathy cluster

# COMMAND ----------

load(file = paste0(raw_data_path, "clusters/risk_test_encephalopathy.rda"))
head(risk_test_encephalopathy)

# COMMAND ----------

# MAGIC %md
# MAGIC ### MCI cluster

# COMMAND ----------

load(file = paste0(raw_data_path, "clusters/risk_test_mci.rda"))
head(risk_test_mci)

# COMMAND ----------

full_risk_result = rbind(risk_test_depression, risk_test_vascular, risk_test_encephalopathy, risk_test_mci) %>% as.data.frame()
filtered_risk_control = full_risk_result %>% 
  # filter(coef > 0 & `p-value` < 0.05) %>%
  mutate(RR = round(as.numeric(`exp(coef)`), 2)) %>% 
  select(start, end, test_group, RR, `p-value`) %>% unique() %>%
  filter(!is.na(RR)) %>%
  filter(test_group %in% c("both vs. either", "both vs. all others", "in order vs. not in order"))
print(dim(filtered_risk_control))
filtered_risk_control

# COMMAND ----------

risk_test_mci %>% filter(test_group == "both vs. either")

# COMMAND ----------

full_risk_result %>% 
  # filter(coef > 0 & `p-value` < 0.05) %>%
  mutate(RR = round(as.numeric(`exp(coef)`), 2)) %>% 
  select(start, end, test_group, RR, `p-value`) %>% unique() %>%
  filter(start == "G31" & end == "I67")

# COMMAND ----------

