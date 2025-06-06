# Databricks notebook source
# MAGIC %md
# MAGIC #### Last updated: 07/02/2024

# COMMAND ----------

library(tidyverse)
'%!in%' = function(x,y)!{'%in%'(x,y)}
raw_data_path = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/mingzhoufu@mednet.ucla.edu/AD_trajectory/data/"

# COMMAND ----------

# MAGIC %md
# MAGIC # 1. Patient demo + EHR characteristics

# COMMAND ----------

load(file = paste0(raw_data_path, "patient_data/mod/AD_demographic_0610.rda"))
print(dim(AD_demographic)) # dim = (24473,15)
# Load in clustered patients
load(file = paste0(raw_data_path, "patient_data/final/pat_traj_final_cluster_0622.rda"))
print(dim(pat_traj_final_cluster)) # dim = (13409,14)

# COMMAND ----------

length(unique(pat_traj_final_cluster$PatientID))

# COMMAND ----------

# MAGIC %md
# MAGIC ### Venn-diagram of clusters

# COMMAND ----------

install.packages("VennDiagram")
library(VennDiagram)

# COMMAND ----------

# patients in multiple clusters
pat_cluster_1 = pat_traj_final_cluster %>% filter(Kmeans == 1) %>% pull(PatientID) %>% unique()
pat_cluster_2 = pat_traj_final_cluster %>% filter(Kmeans == 2) %>% pull(PatientID) %>% unique()
pat_cluster_3 = pat_traj_final_cluster %>% filter(Kmeans == 3) %>% pull(PatientID) %>% unique()
pat_cluster_4 = pat_traj_final_cluster %>% filter(Kmeans == 4) %>% pull(PatientID) %>% unique()

# COMMAND ----------

venn.plot = venn.diagram(
  x = list(
    Set1 = pat_cluster_1,
    Set2 = pat_cluster_2,
    Set3 = pat_cluster_3,
    Set4 = pat_cluster_4
  ),
  category.names = c("Depression", "Vascular", "Encephalopathy", "MCI"),
  filename = NULL,
  output = TRUE,
  fill = c("#26547c", "#ef476f", "#ffd166", "#06d6a0"), # Add colors
  alpha = 0.5, # Set transparency level
  cat.col = c("#26547c", "#ef476f", "#ffd166", "#06d6a0"), # Add colors to category names
  cat.fontface = "bold"
)

# To display the plot
grid.draw(venn.plot)

# COMMAND ----------

# MAGIC %md
# MAGIC ## Unique-cluster patients only

# COMMAND ----------

# overlaps
overlap_pat = c(intersect(pat_cluster_1, pat_cluster_2), intersect(pat_cluster_1, pat_cluster_3), 
              intersect(pat_cluster_1, pat_cluster_4), intersect(pat_cluster_2, pat_cluster_3),
              intersect(pat_cluster_2, pat_cluster_4), intersect(pat_cluster_3, pat_cluster_4)) %>% unique()
print(length(overlap_pat)) # N = 1684

# COMMAND ----------

head(AD_demographic)

# COMMAND ----------

pat_cluster_unique = pat_traj_final_cluster %>% filter(PatientID %!in% overlap_pat)
for (i in 1:4) {
  cluster_i = i
  pat_cluster_i_info = pat_cluster_unique %>% filter(Kmeans == cluster_i) %>% 
    select(PatientID, Kmeans) %>% unique() %>% left_join(AD_demographic) %>%
    mutate(age_AD = as.numeric(AD_Date - BirthDate)/365.25) %>%
    mutate(AD_to_last = as.numeric(date_last_visit - AD_Date)/365.25) %>%
    mutate(AD_to_death = as.numeric(DeathDate - AD_Date)/365.25) %>%
    mutate(first_visit_to_AD = as.numeric(AD_Date - date_first_visit)/365.25) %>%
    mutate(cluster = case_when(
      Kmeans == 1 ~ "Depression",
      Kmeans == 2 ~ "Vascular",
      Kmeans == 3 ~ "Encephalopathy",
      Kmeans == 4 ~ "MCI"
    )) %>% 
    select(PatientID, female, race_ethnicity, PatientLivingStatus, n_icd_3digit, n_encounter, record_length, enc_per_yr, 
          age_last_visit, location_source_value, first_visit_to_AD, age_AD, AD_to_last, AD_to_death, cluster)
  if (i == 1) {
    pat_cluster_unique_info = pat_cluster_i_info
  } else {
    pat_cluster_unique_info = rbind(pat_cluster_unique_info, pat_cluster_i_info)
  }
}
print(dim(pat_cluster_unique_info)) # dim = (4078,15)
pat_cluster_unique_info$cluster = as.factor(pat_cluster_unique_info$cluster)
table(pat_cluster_unique_info$cluster)

# COMMAND ----------

# List of columns to compare
continuous_vars = c("age_last_visit", "record_length", "n_encounter", "n_icd_3digit", "enc_per_yr")
# Function to perform hypothesis tests for continuous variables
continuous_hypothesis_tests = function(data, vars) {
  results = lapply(vars, function(var) {
    wilcox_test = kruskal.test(data[[var]] ~ data[["cluster"]])
    data.frame(
      variable = var,
      p_value = wilcox_test$p.value
    )
  })
  do.call(rbind, results)
}
categorical_vars = c("female", "race_ethnicity", "PatientLivingStatus", "location_source_value")
# Function to summarize categorical variables
# Function to perform hypothesis tests for categorical variables
categorical_hypothesis_tests = function(data, vars) {
  results = lapply(vars, function(var) {
    chisq_test = chisq.test(table(data[[var]], data[["cluster"]]))
    data.frame(
      variable = var,
      p_value = chisq_test$p.value
    )
  })
  do.call(rbind, results)
}

# COMMAND ----------

#======== Continuous variables =============
cont_median = pat_cluster_unique_info %>% group_by(cluster) %>%
    summarise(across(all_of(continuous_vars), list(median = median)), .groups = "drop") %>% t() %>%as.data.frame() 
colnames(cont_median) = cont_median[1, ]
cont_median = cont_median[-1, ]
row.names(cont_median) = NULL
cont_median = cont_median %>% 
  mutate(depression_m = round(as.numeric(Depression), 1), encephalopathy_m = round(as.numeric(Encephalopathy), 1), 
        mci_m = round(as.numeric(MCI), 1), vascular_m = round(as.numeric(Vascular), 1), variable = continuous_vars) %>%
  select(variable, depression_m, encephalopathy_m, mci_m, vascular_m)

cont_IQR = pat_cluster_unique_info %>% group_by(cluster) %>%
    summarise(across(all_of(continuous_vars), list(IQR = IQR)), .groups = "drop") %>% t() %>%as.data.frame() 
colnames(cont_IQR) = cont_IQR[1, ]
cont_IQR = cont_IQR[-1, ]
row.names(cont_IQR) = NULL
cont_IQR = cont_IQR %>% 
  mutate(depression_i = round(as.numeric(Depression), 1), encephalopathy_i = round(as.numeric(Encephalopathy), 1), 
        mci_i = round(as.numeric(MCI), 1), vascular_i = round(as.numeric(Vascular), 1), variable = continuous_vars) %>%
  select(variable, depression_i, encephalopathy_i, mci_i, vascular_i)

cont_tests = continuous_hypothesis_tests(pat_cluster_unique_info, continuous_vars)

cont_full = cont_median %>% full_join(cont_IQR) %>%
  mutate(depression = paste0(depression_m, " (", depression_i, ")"), 
        encephalopathy = paste0(encephalopathy_m, " (", encephalopathy_i, ")"),
        mci = paste0(mci_m, " (", mci_i, ")"), vascular = paste0(vascular_m, " (", vascular_i, ")")) %>%
  select(variable, depression, encephalopathy, mci, vascular) %>% full_join(cont_tests) %>%
  mutate(p_value = as.numeric(p_value)) %>%
  mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value, 3)))) 

# COMMAND ----------

#======== Categorical variables =============
categorical_tests = categorical_hypothesis_tests(pat_cluster_unique_info, categorical_vars)
categorical_tests = categorical_tests %>%
  mutate(p_value = as.numeric(p_value)) %>%
  mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value,3)))) 

# Gender
frequency = table(pat_cluster_unique_info$cluster, pat_cluster_unique_info$female)
formatted_frequency = format(frequency, big.mark = ",")
percentage = prop.table(frequency, margin = 1) * 100
result_table = cbind(formatted_frequency[,2], percentage[,2]) %>% as.data.frame()
names(result_table) = c("n", "percent")
gender_table = result_table %>% mutate(female = paste0(n, " (", round(as.numeric(percent), 1), "%)")) %>%
  select(female) %>% t() %>% as.data.frame() %>% rownames_to_column(var = "variable") %>%
  mutate(p_value = categorical_tests$p_value[1])
names(gender_table) = c("variable", "depression", "encephalopathy", "mci", "vascular", "p_value")

# Race-ethnicity
frequency = table(pat_cluster_unique_info$cluster, pat_cluster_unique_info$race_ethnicity) 
formatted_frequency = format(frequency, big.mark = ",") %>% as.data.frame()
percentage = prop.table(frequency, margin = 1) * 100
formatted_percentage = format(round(percentage, 1)) %>% as.data.frame()
result_table = data.frame()
col_names = names(formatted_percentage)
for (i in 1:(dim(frequency)[2])) {
  if (i == 1) {result_table = cbind(formatted_frequency[,i], formatted_percentage[,i]) %>% as.data.frame()}
  else (result_table = cbind(result_table, formatted_frequency[,i], formatted_percentage[,i]) %>% as.data.frame())
  result_table[[col_names[i]]] = paste0(trimws(result_table[,ncol(result_table)-1]), " (", trimws(result_table[,ncol(result_table)]), "%)")
}
# Select every third column starting from the third column
race_table = result_table[, seq(3, ncol(result_table), by = 3)] %>% t() %>% as.data.frame()
names(race_table) = c("depression", "encephalopathy", "mci", "vascular")
race_table = race_table %>% mutate(p_value = categorical_tests$p_value[2]) %>% rownames_to_column(var = "variable")

# Death
frequency = table(pat_cluster_unique_info$cluster, pat_cluster_unique_info$PatientLivingStatus)
formatted_frequency = format(frequency, big.mark = ",")
percentage = prop.table(frequency, margin = 1) * 100
result_table = cbind(formatted_frequency[,2], percentage[,2]) %>% as.data.frame()
names(result_table) = c("n", "percent")
death_table = result_table %>% mutate(death = paste0(n, " (", round(as.numeric(percent), 1), "%)")) %>%
  select(death) %>% t() %>% as.data.frame() %>% rownames_to_column(var = "variable") %>%
  mutate(p_value = categorical_tests$p_value[3])
names(death_table) = c("variable", "depression", "encephalopathy", "mci", "vascular", "p_value")
# location
frequency = table(pat_cluster_unique_info$cluster, pat_cluster_unique_info$location_source_value) 
formatted_frequency = format(frequency, big.mark = ",") %>% as.data.frame()
percentage = prop.table(frequency, margin = 1) * 100
formatted_percentage = format(round(percentage, 1)) %>% as.data.frame()
result_table = data.frame()
col_names = names(formatted_percentage)
for (i in 1:(dim(frequency)[2])) {
  if (i == 1) {result_table = cbind(formatted_frequency[,i], formatted_percentage[,i]) %>% as.data.frame()}
  else (result_table = cbind(result_table, formatted_frequency[,i], formatted_percentage[,i]) %>% as.data.frame())
  result_table[[col_names[i]]] = paste0(trimws(result_table[,ncol(result_table)-1]), " (", trimws(result_table[,ncol(result_table)]), "%)")
}
# Select every third column starting from the third column
location_table = result_table[, seq(3, ncol(result_table), by = 3)] %>% t() %>% as.data.frame()
names(location_table) = c("depression", "encephalopathy", "mci", "vascular")
location_table = location_table %>% mutate(p_value = categorical_tests$p_value[4]) %>% rownames_to_column(var = "variable")
final_sumstat_unique = rbind(cont_full, gender_table, death_table, c("Race-ethncity", "", "", "", "", ""), 
                race_table, c("Location sources", "", "", "", "", ""), location_table)
final_sumstat_unique

# COMMAND ----------

# List of columns to compare
continuous_vars = c("age_AD", "first_visit_to_AD", "AD_to_last", "AD_to_death")
# Function to perform hypothesis tests for continuous variables
continuous_hypothesis_tests = function(data, vars) {
  results = lapply(vars, function(var) {
    wilcox_test = kruskal.test(data[[var]] ~ data[["cluster"]])
    data.frame(
      variable = var,
      p_value = wilcox_test$p.value
    )
  })
  do.call(rbind, results)
}
#======== Continuous variables =============
cont_median = pat_cluster_unique_info %>% group_by(cluster) %>%
    summarise(across(all_of(continuous_vars), list(median = median), na.rm = T), .groups = "drop") %>% t() %>%as.data.frame() 
colnames(cont_median) = cont_median[1, ]
cont_median = cont_median[-1, ]
row.names(cont_median) = NULL
cont_median = cont_median %>% 
  mutate(depression_m = round(as.numeric(Depression), 1), encephalopathy_m = round(as.numeric(Encephalopathy), 1), 
        mci_m = round(as.numeric(MCI), 1), vascular_m = round(as.numeric(Vascular), 1), variable = continuous_vars) %>%
  select(variable, depression_m, encephalopathy_m, mci_m, vascular_m)

cont_IQR = pat_cluster_unique_info %>% group_by(cluster) %>%
    summarise(across(all_of(continuous_vars), list(IQR = IQR), na.rm = T), .groups = "drop") %>% t() %>%as.data.frame() 
colnames(cont_IQR) = cont_IQR[1, ]
cont_IQR = cont_IQR[-1, ]
row.names(cont_IQR) = NULL
cont_IQR = cont_IQR %>% 
  mutate(depression_i = round(as.numeric(Depression), 1), encephalopathy_i = round(as.numeric(Encephalopathy), 1), 
        mci_i = round(as.numeric(MCI), 1), vascular_i = round(as.numeric(Vascular), 1), variable = continuous_vars) %>%
  select(variable, depression_i, encephalopathy_i, mci_i, vascular_i)

cont_tests = continuous_hypothesis_tests(pat_cluster_unique_info, continuous_vars)
cont_full = cont_median %>% full_join(cont_IQR) %>%
  mutate(depression = paste0(depression_m, " (", depression_i, ")"), 
        encephalopathy = paste0(encephalopathy_m, " (", encephalopathy_i, ")"),
        mci = paste0(mci_m, " (", mci_i, ")"), vascular = paste0(vascular_m, " (", vascular_i, ")")) %>%
  select(variable, depression, encephalopathy, mci, vascular) %>% full_join(cont_tests) %>%
  mutate(p_value = as.numeric(p_value)) %>%
  mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value, 3)))) 
cont_full

# COMMAND ----------

# MAGIC %md
# MAGIC ## Sensitivity: Multi-cluster patients included

# COMMAND ----------

for (i in 1:4) {
  cluster_i = i
  pat_cluster_i_info = pat_traj_final_cluster %>% filter(Kmeans == cluster_i) %>% 
    select(PatientID, Kmeans) %>% unique() %>% left_join(AD_demographic) %>%
    mutate(age_AD = as.numeric(AD_Date - BirthDate)/365.25) %>%
    mutate(AD_to_last = as.numeric(date_last_visit - AD_Date)/365.25) %>%
    mutate(AD_to_death = as.numeric(DeathDate - AD_Date)/365.25) %>%
    mutate(first_visit_to_AD = as.numeric(AD_Date - date_first_visit)/365.25) %>%
    mutate(cluster = case_when(
      Kmeans == 1 ~ "Depression",
      Kmeans == 2 ~ "Vascular",
      Kmeans == 3 ~ "Encephalopathy",
      Kmeans == 4 ~ "MCI"
    )) %>% 
    select(PatientID, female, race_ethnicity, PatientLivingStatus, n_icd_3digit, n_encounter, record_length, enc_per_yr, 
          age_last_visit, location_source_value, first_visit_to_AD, age_AD, AD_to_last, AD_to_death, cluster)
  if (i == 1) {
    pat_cluster_full_info = pat_cluster_i_info
  } else {
    pat_cluster_full_info = rbind(pat_cluster_full_info, pat_cluster_i_info)
  }
}
print(dim(pat_cluster_full_info)) # dim = (7619,15)
pat_cluster_full_info$cluster = as.factor(pat_cluster_full_info$cluster)
table(pat_cluster_full_info$cluster)

# COMMAND ----------

# List of columns to compare
continuous_vars = c("age_last_visit", "record_length", "n_encounter", "n_icd_3digit", "enc_per_yr")
# Function to perform hypothesis tests for continuous variables
continuous_hypothesis_tests = function(data, vars) {
  results = lapply(vars, function(var) {
    wilcox_test = kruskal.test(data[[var]] ~ data[["cluster"]])
    data.frame(
      variable = var,
      p_value = wilcox_test$p.value
    )
  })
  do.call(rbind, results)
}
categorical_vars = c("female", "race_ethnicity", "PatientLivingStatus", "location_source_value")
# Function to summarize categorical variables
# Function to perform hypothesis tests for categorical variables
categorical_hypothesis_tests = function(data, vars) {
  results = lapply(vars, function(var) {
    chisq_test = chisq.test(table(data[[var]], data[["cluster"]]))
    data.frame(
      variable = var,
      p_value = chisq_test$p.value
    )
  })
  do.call(rbind, results)
}

# COMMAND ----------

#======== Continuous variables =============
cont_median = pat_cluster_full_info %>% group_by(cluster) %>%
    summarise(across(all_of(continuous_vars), list(median = median)), .groups = "drop") %>% t() %>%as.data.frame() 
colnames(cont_median) = cont_median[1, ]
cont_median = cont_median[-1, ]
row.names(cont_median) = NULL
cont_median = cont_median %>% 
  mutate(depression_m = round(as.numeric(Depression), 1), encephalopathy_m = round(as.numeric(Encephalopathy), 1), 
        mci_m = round(as.numeric(MCI), 1), vascular_m = round(as.numeric(Vascular), 1), variable = continuous_vars) %>%
  select(variable, depression_m, encephalopathy_m, mci_m, vascular_m)

cont_IQR = pat_cluster_full_info %>% group_by(cluster) %>%
    summarise(across(all_of(continuous_vars), list(IQR = IQR)), .groups = "drop") %>% t() %>%as.data.frame() 
colnames(cont_IQR) = cont_IQR[1, ]
cont_IQR = cont_IQR[-1, ]
row.names(cont_IQR) = NULL
cont_IQR = cont_IQR %>% 
  mutate(depression_i = round(as.numeric(Depression), 1), encephalopathy_i = round(as.numeric(Encephalopathy), 1), 
        mci_i = round(as.numeric(MCI), 1), vascular_i = round(as.numeric(Vascular), 1), variable = continuous_vars) %>%
  select(variable, depression_i, encephalopathy_i, mci_i, vascular_i)

cont_tests = continuous_hypothesis_tests(pat_cluster_unique_info, continuous_vars)

cont_full = cont_median %>% full_join(cont_IQR) %>%
  mutate(depression = paste0(depression_m, " (", depression_i, ")"), 
        encephalopathy = paste0(encephalopathy_m, " (", encephalopathy_i, ")"),
        mci = paste0(mci_m, " (", mci_i, ")"), vascular = paste0(vascular_m, " (", vascular_i, ")")) %>%
  select(variable, depression, encephalopathy, mci, vascular) %>% full_join(cont_tests) %>%
  mutate(p_value = as.numeric(p_value)) %>%
  mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value, 3)))) 

# COMMAND ----------

#======== Categorical variables =============
categorical_tests = categorical_hypothesis_tests(pat_cluster_full_info, categorical_vars)
categorical_tests = categorical_tests %>%
  mutate(p_value = as.numeric(p_value)) %>%
  mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value,3)))) 

# Gender
frequency = table(pat_cluster_full_info$cluster, pat_cluster_full_info$female)
formatted_frequency = format(frequency, big.mark = ",")
percentage = prop.table(frequency, margin = 1) * 100
result_table = cbind(formatted_frequency[,2], percentage[,2]) %>% as.data.frame()
names(result_table) = c("n", "percent")
gender_table = result_table %>% mutate(female = paste0(n, " (", round(as.numeric(percent), 1), "%)")) %>%
  select(female) %>% t() %>% as.data.frame() %>% rownames_to_column(var = "variable") %>%
  mutate(p_value = categorical_tests$p_value[1])
names(gender_table) = c("variable", "depression", "encephalopathy", "mci", "vascular", "p_value")

# Race-ethnicity
frequency = table(pat_cluster_full_info$cluster, pat_cluster_full_info$race_ethnicity) 
formatted_frequency = format(frequency, big.mark = ",") %>% as.data.frame()
percentage = prop.table(frequency, margin = 1) * 100
formatted_percentage = format(round(percentage, 1)) %>% as.data.frame()
result_table = data.frame()
col_names = names(formatted_percentage)
for (i in 1:(dim(frequency)[2])) {
  if (i == 1) {result_table = cbind(formatted_frequency[,i], formatted_percentage[,i]) %>% as.data.frame()}
  else (result_table = cbind(result_table, formatted_frequency[,i], formatted_percentage[,i]) %>% as.data.frame())
  result_table[[col_names[i]]] = paste0(trimws(result_table[,ncol(result_table)-1]), " (", trimws(result_table[,ncol(result_table)]), "%)")
}
# Select every third column starting from the third column
race_table = result_table[, seq(3, ncol(result_table), by = 3)] %>% t() %>% as.data.frame()
names(race_table) = c("depression", "encephalopathy", "mci", "vascular")
race_table = race_table %>% mutate(p_value = categorical_tests$p_value[2]) %>% rownames_to_column(var = "variable")

# Death
frequency = table(pat_cluster_full_info$cluster, pat_cluster_full_info$PatientLivingStatus)
formatted_frequency = format(frequency, big.mark = ",")
percentage = prop.table(frequency, margin = 1) * 100
result_table = cbind(formatted_frequency[,2], percentage[,2]) %>% as.data.frame()
names(result_table) = c("n", "percent")
death_table = result_table %>% mutate(death = paste0(n, " (", round(as.numeric(percent), 1), "%)")) %>%
  select(death) %>% t() %>% as.data.frame() %>% rownames_to_column(var = "variable") %>%
  mutate(p_value = categorical_tests$p_value[3])
names(death_table) = c("variable", "depression", "encephalopathy", "mci", "vascular", "p_value")
# location
frequency = table(pat_cluster_full_info$cluster, pat_cluster_full_info$location_source_value) 
formatted_frequency = format(frequency, big.mark = ",") %>% as.data.frame()
percentage = prop.table(frequency, margin = 1) * 100
formatted_percentage = format(round(percentage, 1)) %>% as.data.frame()
result_table = data.frame()
col_names = names(formatted_percentage)
for (i in 1:(dim(frequency)[2])) {
  if (i == 1) {result_table = cbind(formatted_frequency[,i], formatted_percentage[,i]) %>% as.data.frame()}
  else (result_table = cbind(result_table, formatted_frequency[,i], formatted_percentage[,i]) %>% as.data.frame())
  result_table[[col_names[i]]] = paste0(trimws(result_table[,ncol(result_table)-1]), " (", trimws(result_table[,ncol(result_table)]), "%)")
}
# Select every third column starting from the third column
location_table = result_table[, seq(3, ncol(result_table), by = 3)] %>% t() %>% as.data.frame()
names(location_table) = c("depression", "encephalopathy", "mci", "vascular")
location_table = location_table %>% mutate(p_value = categorical_tests$p_value[4]) %>% rownames_to_column(var = "variable")
final_sumstat_full = rbind(cont_full, gender_table, death_table, c("Race-ethncity", "", "", "", "", ""), 
                race_table, c("Location sources", "", "", "", "", ""), location_table)
final_sumstat_full

# COMMAND ----------

# List of columns to compare
continuous_vars = c("age_AD", "first_visit_to_AD", "AD_to_last", "AD_to_death")
# Function to perform hypothesis tests for continuous variables
continuous_hypothesis_tests = function(data, vars) {
  results = lapply(vars, function(var) {
    wilcox_test = kruskal.test(data[[var]] ~ data[["cluster"]])
    data.frame(
      variable = var,
      p_value = wilcox_test$p.value
    )
  })
  do.call(rbind, results)
}
#======== Continuous variables =============
cont_median = pat_cluster_full_info %>% group_by(cluster) %>%
    summarise(across(all_of(continuous_vars), list(median = median), na.rm = T), .groups = "drop") %>% t() %>%as.data.frame() 
colnames(cont_median) = cont_median[1, ]
cont_median = cont_median[-1, ]
row.names(cont_median) = NULL
cont_median = cont_median %>% 
  mutate(depression_m = round(as.numeric(Depression), 1), encephalopathy_m = round(as.numeric(Encephalopathy), 1), 
        mci_m = round(as.numeric(MCI), 1), vascular_m = round(as.numeric(Vascular), 1), variable = continuous_vars) %>%
  select(variable, depression_m, encephalopathy_m, mci_m, vascular_m)

cont_IQR = pat_cluster_full_info %>% group_by(cluster) %>%
    summarise(across(all_of(continuous_vars), list(IQR = IQR), na.rm = T), .groups = "drop") %>% t() %>%as.data.frame() 
colnames(cont_IQR) = cont_IQR[1, ]
cont_IQR = cont_IQR[-1, ]
row.names(cont_IQR) = NULL
cont_IQR = cont_IQR %>% 
  mutate(depression_i = round(as.numeric(Depression), 1), encephalopathy_i = round(as.numeric(Encephalopathy), 1), 
        mci_i = round(as.numeric(MCI), 1), vascular_i = round(as.numeric(Vascular), 1), variable = continuous_vars) %>%
  select(variable, depression_i, encephalopathy_i, mci_i, vascular_i)

cont_tests = continuous_hypothesis_tests(pat_cluster_full_info, continuous_vars)
cont_full = cont_median %>% full_join(cont_IQR) %>%
  mutate(depression = paste0(depression_m, " (", depression_i, ")"), 
        encephalopathy = paste0(encephalopathy_m, " (", encephalopathy_i, ")"),
        mci = paste0(mci_m, " (", mci_i, ")"), vascular = paste0(vascular_m, " (", vascular_i, ")")) %>%
  select(variable, depression, encephalopathy, mci, vascular) %>% full_join(cont_tests) %>%
  mutate(p_value = as.numeric(p_value)) %>%
  mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value, 3)))) 
cont_full

# COMMAND ----------

# MAGIC %md
# MAGIC # 2. Symptoms data

# COMMAND ----------

load(file = paste0(raw_data_path, "patient_data/mod/AD_encounters_filter_0610.rda"))
print(dim(AD_encounters_filter)) # dim = (974718,5)

# COMMAND ----------

symptom_encounters = AD_encounters_filter %>% mutate(ICD_3digit = substr(ICD, 1, 3)) %>%
  select(-c(ICD, date_first_visit)) %>% unique() %>% group_by(PatientID, ICD_3digit) %>%
  arrange(EncounterDate) %>% slice(1) %>% ungroup() %>%
  mutate(chapter = substr(ICD_3digit, 1, 1)) %>%
  filter(chapter == "R") %>% select(PatientID, EncounterDate, ICD_3digit, AD_Date) %>% unique()
print(dim(symptom_encounters)) # dim = (128393,4)

# COMMAND ----------

# Function to perform hypothesis tests for categorical variables
categorical_hypothesis_tests = function(data, vars) {
  results = lapply(vars, function(var) {
    chisq_test = chisq.test(table(data[[var]], data[["cluster"]]))
    data.frame(
      variable = var,
      p_value = chisq_test$p.value
    )
  })
  do.call(rbind, results)
}
# Function to count number of rows with values > 0 for each column
count_greater_than_zero = function(df) {
  sapply(df, function(column) sum(column > 0))
}

# COMMAND ----------

# MAGIC %md
# MAGIC ### Cumulated symptoms (of all time)

# COMMAND ----------

# MAGIC %md
# MAGIC #### Cognition symptoms

# COMMAND ----------

symptom_cog = symptom_encounters %>% filter(EncounterDate <= AD_Date) %>%
  mutate(subchapter = substr(ICD_3digit, 1, 2)) %>% filter(subchapter == "R4") %>% 
  select(PatientID, ICD_3digit) %>% mutate(value = 1) %>%
  pivot_wider(names_from = ICD_3digit, values_from = value, values_fill = 0)
cluster_symptom_cog = pat_cluster_unique_info %>% select(PatientID, cluster) %>% unique() %>% left_join(symptom_cog)
# Fill NAs with 0
cluster_symptom_cog[is.na(cluster_symptom_cog)] = 0
print(dim(cluster_symptom_cog)) # dim = (4078,12)
# Calculate column sums excluding the first column
column_sums = colSums(cluster_symptom_cog[, -(1:2)])
# Print the column sums
print(sort(column_sums))

# COMMAND ----------

symptom_cog_freq = column_sums %>% as.data.frame() %>% rownames_to_column()
names(symptom_cog_freq) = c("ICD", "freq")
symptom_cog_freq = symptom_cog_freq %>% arrange(ICD)
symptom_of_interest = setdiff(names(cluster_symptom_cog), c("PatientID", "cluster")) %>% unique() %>% sort()
symptom_subset = cluster_symptom_cog %>% select(PatientID, cluster, all_of(symptom_of_interest)) 
categorical_tests = categorical_hypothesis_tests(symptom_subset, symptom_of_interest) %>% mutate(p_value = as.numeric(p_value)) 
symptom_cog_sig = symptom_cog_freq %>% inner_join(categorical_tests, by = c("ICD" = "variable")) %>%
  mutate(fdr_adj_p = p.adjust(p_value, method = "fdr")) 
print(dim(symptom_cog_sig)) # dim = c(10,4)
head(symptom_cog_sig, 10)

# COMMAND ----------

symptom_of_interest = c("R40", "R41", "R42", "R43", "R44", "R45", "R46")
for (i in 1:length(symptom_of_interest)) {
  symptom = symptom_of_interest[i]
  frequency = table(symptom_subset$cluster, symptom_subset[[symptom]])
  formatted_frequency = format(frequency, big.mark = ",")
  percentage = prop.table(frequency, margin = 1) * 100
  result_table = cbind(formatted_frequency[,2], percentage[,2]) %>% as.data.frame()
  names(result_table) = c("n", "percent")
  final_table = result_table %>% mutate(test_symptom = paste0(n, " (", round(as.numeric(percent), 1), "%)")) %>%
    select(test_symptom) %>% t() %>% as.data.frame() %>% rownames_to_column(var = "variable") %>%
    mutate(p_value = symptom_cog_sig$fdr_adj_p[i]) %>% 
    mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value,3)))) 
  names(final_table) = c("variable", "depression", "encephalopathy", "mci", "vascular", "adj_p_value")
  final_table$variable = symptom
  if (i == 1) {
    symptom_summary_tbl = final_table
  } else {
    symptom_summary_tbl = rbind(symptom_summary_tbl, final_table)
  }
}
symptom_summary_tbl

# COMMAND ----------

# MAGIC %md
# MAGIC #### Other symptoms (subchapter)

# COMMAND ----------

symptom_other = symptom_encounters %>% filter(EncounterDate <= AD_Date) %>%
  mutate(subchapter = case_when(
    substr(ICD_3digit, 1, 2) == "R0" ~ "circulatory",
    substr(ICD_3digit, 1, 2) == "R1" ~ "digestive",
    ICD_3digit %in% c("R20", "R21", "R22", "R23") ~ "skin",
    ICD_3digit %in% c("R25", "R26", "R27", "R28", "R29") ~ "nerve_muscul",
    substr(ICD_3digit, 1, 2) == "R3" ~ "urinary",
    ICD_3digit %in% c("R40", "R41", "R42", "R43", "R44", "R45", "R46") ~ "cognition",
    ICD_3digit %in% c("R47", "R48", "R49") ~ "speech",
    substr(ICD_3digit, 1, 2) %in% c("R5", "R6") ~ "general",
    substr(ICD_3digit, 1, 2) == "R7" ~ "exam_blood",
    ICD_3digit %in% c("R80", "R81", "R82") ~ "exam_urine",
    ICD_3digit %in% c("R83", "R84", "R85", "R86", "R87", "R88", "R89") ~ "exam_other",
    ICD_3digit %in% c("R90", "R91", "R92", "R93", "R94") ~ "imaging",
    ICD_3digit %in% c("R95", "R96", "R97", "R98", "R99") ~ "mortality"
  )) %>% group_by(PatientID, subchapter) %>% mutate(n_subchapter = n_distinct(ICD_3digit)) %>% ungroup() %>%
  select(PatientID, subchapter, n_subchapter) %>% unique() %>%
  pivot_wider(names_from = subchapter, values_from = n_subchapter, values_fill = 0)
cluster_symptom_other = pat_cluster_unique_info %>% select(PatientID, cluster) %>% unique() %>% left_join(symptom_other)
# Fill NAs with 0
cluster_symptom_other[is.na(cluster_symptom_other)] = 0
print(dim(cluster_symptom_other)) # dim = (4078,15)
# Calculate column sums excluding the first column
column_sums = colSums(cluster_symptom_other[, -(1:2)])
# Print the column sums
print(sort(column_sums))
non_zero_pat = count_greater_than_zero(cluster_symptom_other[, -(1:2)]) %>% as.data.frame() %>% rownames_to_column()
names(non_zero_pat) = c("ICD", "n_pat")

# COMMAND ----------

symptom_other_freq = column_sums %>% as.data.frame() %>% rownames_to_column()
names(symptom_other_freq) = c("ICD", "freq")
symptom_other_freq = symptom_other_freq %>% arrange(ICD)
symptom_of_interest = setdiff(names(cluster_symptom_other), c("PatientID", "cluster")) %>% unique() %>% sort()
symptom_subset = cluster_symptom_other %>% select(PatientID, cluster, all_of(symptom_of_interest)) 
categorical_tests = categorical_hypothesis_tests(symptom_subset, symptom_of_interest) %>% mutate(p_value = as.numeric(p_value)) 
symptom_other_sig = symptom_other_freq %>% inner_join(categorical_tests, by = c("ICD" = "variable")) %>%
  mutate(fdr_adj_p = p.adjust(p_value, method = "fdr")) %>% left_join(non_zero_pat) %>% 
  mutate(n_per_pat = freq/n_pat) %>% select(ICD, freq, n_pat, n_per_pat, fdr_adj_p)
print(dim(symptom_other_sig)) # dim = c(13,4)
head(symptom_other_sig, 10)

# COMMAND ----------

symptom_of_interest = symptom_other_sig$ICD
for (i in 1:length(symptom_of_interest)) {
  symptom = symptom_of_interest[i]
  frequency = table(symptom_subset$cluster, symptom_subset[[symptom]])
  formatted_frequency = format(frequency, big.mark = ",")
  percentage = prop.table(frequency, margin = 1) * 100
  result_table = cbind(formatted_frequency[,2], percentage[,2]) %>% as.data.frame()
  names(result_table) = c("n", "percent")
  final_table = result_table %>% mutate(test_symptom = paste0(n, " (", round(as.numeric(percent), 1), "%)")) %>%
    select(test_symptom) %>% t() %>% as.data.frame() %>% rownames_to_column(var = "variable") %>%
    mutate(p_value = symptom_other_sig$fdr_adj_p[i]) %>% 
    mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value,3)))) 
  names(final_table) = c("variable", "depression", "encephalopathy", "mci", "vascular", "adj_p_value")
  final_table$variable = symptom
  if (i == 1) {
    symptom_summary_tbl = final_table
  } else {
    symptom_summary_tbl = rbind(symptom_summary_tbl, final_table)
  }
}
symptom_summary_tbl

# COMMAND ----------

# MAGIC %md
# MAGIC ### 3-years before AD diagnosis

# COMMAND ----------

# MAGIC %md
# MAGIC #### Cognition symptoms

# COMMAND ----------

symptom_cog_3yr = symptom_encounters %>% filter(EncounterDate >= (AD_Date - lubridate::years(3))) %>%
  mutate(subchapter = substr(ICD_3digit, 1, 2)) %>% filter(subchapter == "R4") %>% 
  select(PatientID, ICD_3digit) %>% mutate(value = 1) %>%
  pivot_wider(names_from = ICD_3digit, values_from = value, values_fill = 0)
cluster_symptom_cog = pat_cluster_unique_info %>% select(PatientID, cluster) %>% unique() %>% left_join(symptom_cog_3yr)
# Fill NAs with 0
cluster_symptom_cog[is.na(cluster_symptom_cog)] = 0
print(dim(cluster_symptom_cog)) # dim = (4078,12)
# Calculate column sums excluding the first column
column_sums = colSums(cluster_symptom_cog[, -(1:2)])
# Print the column sums
print(sort(column_sums))

# COMMAND ----------

symptom_cog_freq = column_sums %>% as.data.frame() %>% rownames_to_column()
names(symptom_cog_freq) = c("ICD", "freq")
symptom_cog_freq = symptom_cog_freq %>% arrange(ICD)
symptom_of_interest = setdiff(names(cluster_symptom_cog), c("PatientID", "cluster")) %>% unique() %>% sort()
symptom_subset = cluster_symptom_cog %>% select(PatientID, cluster, all_of(symptom_of_interest)) 
categorical_tests = categorical_hypothesis_tests(symptom_subset, symptom_of_interest) %>% mutate(p_value = as.numeric(p_value)) 
symptom_cog_sig = symptom_cog_freq %>% inner_join(categorical_tests, by = c("ICD" = "variable")) %>%
  mutate(fdr_adj_p = p.adjust(p_value, method = "fdr")) 
print(dim(symptom_cog_sig)) # dim = c(10,4)
head(symptom_cog_sig, 10)

# COMMAND ----------

symptom_of_interest = c("R40", "R41", "R42", "R43", "R44", "R45", "R46")
for (i in 1:length(symptom_of_interest)) {
  symptom = symptom_of_interest[i]
  frequency = table(symptom_subset$cluster, symptom_subset[[symptom]])
  formatted_frequency = format(frequency, big.mark = ",")
  percentage = prop.table(frequency, margin = 1) * 100
  result_table = cbind(formatted_frequency[,2], percentage[,2]) %>% as.data.frame()
  names(result_table) = c("n", "percent")
  final_table = result_table %>% mutate(test_symptom = paste0(n, " (", round(as.numeric(percent), 1), "%)")) %>%
    select(test_symptom) %>% t() %>% as.data.frame() %>% rownames_to_column(var = "variable") %>%
    mutate(p_value = symptom_cog_sig$fdr_adj_p[i]) %>% 
    mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value,3)))) 
  names(final_table) = c("variable", "depression", "encephalopathy", "mci", "vascular", "adj_p_value")
  final_table$variable = symptom
  if (i == 1) {
    symptom_summary_tbl = final_table
  } else {
    symptom_summary_tbl = rbind(symptom_summary_tbl, final_table)
  }
}
symptom_summary_tbl

# COMMAND ----------

# MAGIC %md
# MAGIC #### Other symptoms (subchapter)

# COMMAND ----------

symptom_other_3yr = symptom_encounters %>% filter(EncounterDate >= (AD_Date - lubridate::years(3))) %>%
  mutate(subchapter = case_when(
    substr(ICD_3digit, 1, 2) == "R0" ~ "circulatory",
    substr(ICD_3digit, 1, 2) == "R1" ~ "digestive",
    ICD_3digit %in% c("R20", "R21", "R22", "R23") ~ "skin",
    ICD_3digit %in% c("R25", "R26", "R27", "R28", "R29") ~ "nerve_muscul",
    substr(ICD_3digit, 1, 2) == "R3" ~ "urinary",
    ICD_3digit %in% c("R40", "R41", "R42", "R43", "R44", "R45", "R46") ~ "cognition",
    ICD_3digit %in% c("R47", "R48", "R49") ~ "speech",
    substr(ICD_3digit, 1, 2) %in% c("R5", "R6") ~ "general",
    substr(ICD_3digit, 1, 2) == "R7" ~ "exam_blood",
    ICD_3digit %in% c("R80", "R81", "R82") ~ "exam_urine",
    ICD_3digit %in% c("R83", "R84", "R85", "R86", "R87", "R88", "R89") ~ "exam_other",
    ICD_3digit %in% c("R90", "R91", "R92", "R93", "R94") ~ "imaging",
    ICD_3digit %in% c("R95", "R96", "R97", "R98", "R99") ~ "mortality"
  )) %>% group_by(PatientID, subchapter) %>% mutate(n_subchapter = n_distinct(ICD_3digit)) %>% ungroup() %>%
  select(PatientID, subchapter, n_subchapter) %>% unique() %>%
  pivot_wider(names_from = subchapter, values_from = n_subchapter, values_fill = 0)
cluster_symptom_other = pat_cluster_unique_info %>% select(PatientID, cluster) %>% unique() %>% left_join(symptom_other_3yr)
# Fill NAs with 0
cluster_symptom_other[is.na(cluster_symptom_other)] = 0
print(dim(cluster_symptom_other)) # dim = (4078,15)
# Calculate column sums excluding the first column
column_sums = colSums(cluster_symptom_other[, -(1:2)])
# Print the column sums
print(sort(column_sums))
non_zero_pat = count_greater_than_zero(cluster_symptom_other[, -(1:2)]) %>% as.data.frame() %>% rownames_to_column()
names(non_zero_pat) = c("ICD", "n_pat")

# COMMAND ----------

symptom_other_freq = column_sums %>% as.data.frame() %>% rownames_to_column()
names(symptom_other_freq) = c("ICD", "freq")
symptom_other_freq = symptom_other_freq %>% arrange(ICD)
symptom_of_interest = setdiff(names(cluster_symptom_other), c("PatientID", "cluster")) %>% unique() %>% sort()
symptom_subset = cluster_symptom_other %>% select(PatientID, cluster, all_of(symptom_of_interest)) 
categorical_tests = categorical_hypothesis_tests(symptom_subset, symptom_of_interest) %>% mutate(p_value = as.numeric(p_value)) 
symptom_other_sig = symptom_other_freq %>% inner_join(categorical_tests, by = c("ICD" = "variable")) %>%
  mutate(fdr_adj_p = p.adjust(p_value, method = "fdr")) %>% left_join(non_zero_pat) %>% 
  mutate(n_per_pat = freq/n_pat) %>% select(ICD, freq, n_pat, n_per_pat, fdr_adj_p)
print(dim(symptom_other_sig)) # dim = c(13,4)
head(symptom_other_sig, 10)

# COMMAND ----------

symptom_of_interest = symptom_other_sig$ICD
for (i in 1:length(symptom_of_interest)) {
  symptom = symptom_of_interest[i]
  frequency = table(symptom_subset$cluster, symptom_subset[[symptom]])
  formatted_frequency = format(frequency, big.mark = ",")
  percentage = prop.table(frequency, margin = 1) * 100
  result_table = cbind(formatted_frequency[,2], percentage[,2]) %>% as.data.frame()
  names(result_table) = c("n", "percent")
  final_table = result_table %>% mutate(test_symptom = paste0(n, " (", round(as.numeric(percent), 1), "%)")) %>%
    select(test_symptom) %>% t() %>% as.data.frame() %>% rownames_to_column(var = "variable") %>%
    mutate(p_value = symptom_other_sig$fdr_adj_p[i]) %>% 
    mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value,3)))) 
  names(final_table) = c("variable", "depression", "encephalopathy", "mci", "vascular", "adj_p_value")
  final_table$variable = symptom
  if (i == 1) {
    symptom_summary_tbl = final_table
  } else {
    symptom_summary_tbl = rbind(symptom_summary_tbl, final_table)
  }
}
symptom_summary_tbl

# COMMAND ----------

# MAGIC %md
# MAGIC ### 1-year before AD diagnosis

# COMMAND ----------

# MAGIC %md
# MAGIC #### Cognition symptoms

# COMMAND ----------

symptom_cog_1yr = symptom_encounters %>% filter(EncounterDate >= (AD_Date - lubridate::years(1))) %>%
  mutate(subchapter = substr(ICD_3digit, 1, 2)) %>% filter(subchapter == "R4") %>% 
  select(PatientID, ICD_3digit) %>% mutate(value = 1) %>%
  pivot_wider(names_from = ICD_3digit, values_from = value, values_fill = 0)
cluster_symptom_cog = pat_cluster_unique_info %>% select(PatientID, cluster) %>% unique() %>% left_join(symptom_cog_1yr)
# Fill NAs with 0
cluster_symptom_cog[is.na(cluster_symptom_cog)] = 0
print(dim(cluster_symptom_cog)) # dim = (4078,12)
# Calculate column sums excluding the first column
column_sums = colSums(cluster_symptom_cog[, -(1:2)])
# Print the column sums
print(sort(column_sums))

# COMMAND ----------

symptom_cog_freq = column_sums %>% as.data.frame() %>% rownames_to_column()
names(symptom_cog_freq) = c("ICD", "freq")
symptom_cog_freq = symptom_cog_freq %>% arrange(ICD)
symptom_of_interest = setdiff(names(cluster_symptom_cog), c("PatientID", "cluster")) %>% unique() %>% sort()
symptom_subset = cluster_symptom_cog %>% select(PatientID, cluster, all_of(symptom_of_interest)) 
categorical_tests = categorical_hypothesis_tests(symptom_subset, symptom_of_interest) %>% mutate(p_value = as.numeric(p_value)) 
symptom_cog_sig = symptom_cog_freq %>% inner_join(categorical_tests, by = c("ICD" = "variable")) %>%
  mutate(fdr_adj_p = p.adjust(p_value, method = "fdr")) 
print(dim(symptom_cog_sig)) # dim = c(10,4)
head(symptom_cog_sig, 10)

# COMMAND ----------

symptom_of_interest = c("R40", "R41", "R42", "R43", "R44", "R45", "R46")
for (i in 1:length(symptom_of_interest)) {
  symptom = symptom_of_interest[i]
  frequency = table(symptom_subset$cluster, symptom_subset[[symptom]])
  formatted_frequency = format(frequency, big.mark = ",")
  percentage = prop.table(frequency, margin = 1) * 100
  result_table = cbind(formatted_frequency[,2], percentage[,2]) %>% as.data.frame()
  names(result_table) = c("n", "percent")
  final_table = result_table %>% mutate(test_symptom = paste0(n, " (", round(as.numeric(percent), 1), "%)")) %>%
    select(test_symptom) %>% t() %>% as.data.frame() %>% rownames_to_column(var = "variable") %>%
    mutate(p_value = symptom_cog_sig$fdr_adj_p[i]) %>% 
    mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value,3)))) 
  names(final_table) = c("variable", "depression", "encephalopathy", "mci", "vascular", "adj_p_value")
  final_table$variable = symptom
  if (i == 1) {
    symptom_summary_tbl = final_table
  } else {
    symptom_summary_tbl = rbind(symptom_summary_tbl, final_table)
  }
}
symptom_summary_tbl

# COMMAND ----------

# MAGIC %md
# MAGIC #### Other symptoms (subchapter)

# COMMAND ----------

symptom_other_1yr = symptom_encounters %>% filter(EncounterDate >= (AD_Date - lubridate::years(1))) %>%
  mutate(subchapter = case_when(
    substr(ICD_3digit, 1, 2) == "R0" ~ "circulatory",
    substr(ICD_3digit, 1, 2) == "R1" ~ "digestive",
    ICD_3digit %in% c("R20", "R21", "R22", "R23") ~ "skin",
    ICD_3digit %in% c("R25", "R26", "R27", "R28", "R29") ~ "nerve_muscul",
    substr(ICD_3digit, 1, 2) == "R3" ~ "urinary",
    ICD_3digit %in% c("R40", "R41", "R42", "R43", "R44", "R45", "R46") ~ "cognition",
    ICD_3digit %in% c("R47", "R48", "R49") ~ "speech",
    substr(ICD_3digit, 1, 2) %in% c("R5", "R6") ~ "general",
    substr(ICD_3digit, 1, 2) == "R7" ~ "exam_blood",
    ICD_3digit %in% c("R80", "R81", "R82") ~ "exam_urine",
    ICD_3digit %in% c("R83", "R84", "R85", "R86", "R87", "R88", "R89") ~ "exam_other",
    ICD_3digit %in% c("R90", "R91", "R92", "R93", "R94") ~ "imaging",
    ICD_3digit %in% c("R95", "R96", "R97", "R98", "R99") ~ "mortality"
  )) %>% group_by(PatientID, subchapter) %>% mutate(n_subchapter = n_distinct(ICD_3digit)) %>% ungroup() %>%
  select(PatientID, subchapter, n_subchapter) %>% unique() %>%
  pivot_wider(names_from = subchapter, values_from = n_subchapter, values_fill = 0)
cluster_symptom_other = pat_cluster_unique_info %>% select(PatientID, cluster) %>% unique() %>% left_join(symptom_other_1yr)
# Fill NAs with 0
cluster_symptom_other[is.na(cluster_symptom_other)] = 0
print(dim(cluster_symptom_other)) # dim = (4078,15)
# Calculate column sums excluding the first column
column_sums = colSums(cluster_symptom_other[, -(1:2)])
# Print the column sums
print(sort(column_sums))
non_zero_pat = count_greater_than_zero(cluster_symptom_other[, -(1:2)]) %>% as.data.frame() %>% rownames_to_column()
names(non_zero_pat) = c("ICD", "n_pat")

# COMMAND ----------

symptom_other_freq = column_sums %>% as.data.frame() %>% rownames_to_column()
names(symptom_other_freq) = c("ICD", "freq")
symptom_other_freq = symptom_other_freq %>% arrange(ICD)
symptom_of_interest = setdiff(names(cluster_symptom_other), c("PatientID", "cluster")) %>% unique() %>% sort()
symptom_subset = cluster_symptom_other %>% select(PatientID, cluster, all_of(symptom_of_interest)) 
categorical_tests = categorical_hypothesis_tests(symptom_subset, symptom_of_interest) %>% mutate(p_value = as.numeric(p_value)) 
symptom_other_sig = symptom_other_freq %>% inner_join(categorical_tests, by = c("ICD" = "variable")) %>%
  mutate(fdr_adj_p = p.adjust(p_value, method = "fdr")) %>% left_join(non_zero_pat) %>% 
  mutate(n_per_pat = freq/n_pat) %>% select(ICD, freq, n_pat, n_per_pat, fdr_adj_p)
print(dim(symptom_other_sig)) # dim = c(13,4)
head(symptom_other_sig, 10)

# COMMAND ----------

symptom_of_interest = symptom_other_sig$ICD
for (i in 1:length(symptom_of_interest)) {
  symptom = symptom_of_interest[i]
  frequency = table(symptom_subset$cluster, symptom_subset[[symptom]])
  formatted_frequency = format(frequency, big.mark = ",")
  percentage = prop.table(frequency, margin = 1) * 100
  result_table = cbind(formatted_frequency[,2], percentage[,2]) %>% as.data.frame()
  names(result_table) = c("n", "percent")
  final_table = result_table %>% mutate(test_symptom = paste0(n, " (", round(as.numeric(percent), 1), "%)")) %>%
    select(test_symptom) %>% t() %>% as.data.frame() %>% rownames_to_column(var = "variable") %>%
    mutate(p_value = symptom_other_sig$fdr_adj_p[i]) %>% 
    mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value,3)))) 
  names(final_table) = c("variable", "depression", "encephalopathy", "mci", "vascular", "adj_p_value")
  final_table$variable = symptom
  if (i == 1) {
    symptom_summary_tbl = final_table
  } else {
    symptom_summary_tbl = rbind(symptom_summary_tbl, final_table)
  }
}
symptom_summary_tbl

# COMMAND ----------

# MAGIC %md
# MAGIC # 3. Time to AD and time to death

# COMMAND ----------

# MAGIC %md
# MAGIC ### Kaplan-Meier survival curve

# COMMAND ----------

install.packages("survival")
library(remotes)
install_version("Matrix", version = "1.6-5")
install.packages("survminer")
library(survival)
library(ggplot2)
library(survminer)

# COMMAND ----------

names(AD_demographic)

# COMMAND ----------

length(unique(pat_cluster_unique$PatientID))

# COMMAND ----------

# MAGIC %md
# MAGIC ### Calculate first date in the trajectory 

# COMMAND ----------

load(file = paste0(raw_data_path, "FG_test/other_risk/diag_beforeAD_filter_0610.rda"))
print(dim(diag_beforeAD_filter)) # dim = (393601,4)
head(diag_beforeAD_filter)

# COMMAND ----------

for (i in 1:4) {
  cluster_i = i
  pat_cluster_id = pat_cluster_unique %>% filter(Kmeans == cluster_i) %>% 
    pull(PatientID) %>% unique()
  # print(length(pat_cluster_id))
  load(file = paste0(raw_data_path, "clusters/filter_edges_", cluster_i, ".rda"))
  nodes_in_network_i = unique(c(filter_edges$start, filter_edges$end))
  pat_nodes_in_network = diag_beforeAD_filter %>% filter(PatientID %in% pat_cluster_id) %>%
    filter(ICD_3digit %in% nodes_in_network_i) %>% group_by(PatientID) %>%
    summarize(nodeDate = min(EncounterDate))
  if (i == 1) {
    pat_node_info = pat_nodes_in_network
  } else {
    pat_node_info = rbind(pat_node_info, pat_nodes_in_network)
  }
}
print(dim(pat_node_info))
head(pat_node_info)

# COMMAND ----------

nodes_in_network_i

# COMMAND ----------

# pat_cluster_unique = pat_traj_final_cluster %>% filter(PatientID %!in% overlap_pat)
for (i in 1:4) {
  cluster_i = i
  pat_cluster_i_info = pat_cluster_unique %>% filter(Kmeans == cluster_i) %>% 
    select(PatientID, Kmeans) %>% unique() %>% left_join(AD_demographic) %>%
    mutate(age_AD = as.numeric(AD_Date - BirthDate)/365.25,
          if_death = if_else(PatientLivingStatus == "Deceased", 1, 0),
          age_first_visit = age_last_visit - record_length,
          first_visit_to_AD = as.numeric(AD_Date - date_first_visit)/365.25,
          AD_to_last_visit = as.numeric(date_last_visit - AD_Date)/365.25,
          AD_to_death = as.numeric(DeathDate - AD_Date)/365.25) %>%
    mutate(cluster = case_when(
      Kmeans == 1 ~ "Depression",
      Kmeans == 2 ~ "Vascular",
      Kmeans == 3 ~ "Encephalopathy",
      Kmeans == 4 ~ "MCI"
    )) %>% left_join(pat_node_info) %>%
    mutate(first_node_to_AD = as.numeric(AD_Date - nodeDate)/365.25) %>%
    select(PatientID, AD_Date, age_AD, first_visit_to_AD, first_node_to_AD, AD_to_last_visit, if_death, 
    AD_to_death, age_first_visit, age_last_visit, cluster)
  if (i == 1) {
    pat_cluster_unique_info = pat_cluster_i_info
  } else {
    pat_cluster_unique_info = rbind(pat_cluster_unique_info, pat_cluster_i_info)
  }
}
print(dim(pat_cluster_unique_info)) # dim = (4078,11)
pat_cluster_unique_info$cluster = as.factor(pat_cluster_unique_info$cluster)
table(pat_cluster_unique_info$cluster, pat_cluster_unique_info$if_death)

# COMMAND ----------

# MAGIC %md
# MAGIC #### Time to AD (event)

# COMMAND ----------

time_to_AD = pat_cluster_unique_info %>% select(PatientID, first_node_to_AD, cluster) %>%
  mutate(status = 1) %>% rename(time = first_node_to_AD) 
# Fit Kaplan-Meier curves by groups
km_fit = survfit(Surv(time, status) ~ cluster, data = time_to_AD)

# COMMAND ----------

# Plot Kaplan-Meier curves
ggsurvplot(km_fit, conf.int = TRUE, pval = TRUE, risk.table = TRUE, censor.size = 0.6, size = 0.5,
          legend.title = "Cluster", legend.labs = c("Depression", "Encephalopathy", "MCI", "Vascular"),
          palette = c("#26547c", "#ffd166", "#06d6a0", "#ef476f"), 
          title = "Kaplan-Meier Curve for first node to AD Survival", 
          risk.table.height = .3)

# COMMAND ----------

# MAGIC %md
# MAGIC #### AD to death

# COMMAND ----------

AD_to_death = pat_cluster_unique_info %>% 
  mutate(time = case_when(
    if_death == 1 ~ AD_to_death,
    if_death == 0 ~ AD_to_last_visit
  )) %>% mutate(time = as.numeric(time)) %>%
  select(PatientID, time, cluster, if_death) %>% 
  mutate(status = case_when(
    if_death == 1 ~ 1,
    if_death == 0 ~ 0
  )) 
# Fit Kaplan-Meier curves by groups
km_fit = survfit(Surv(time, status) ~ cluster, data = AD_to_death)

# COMMAND ----------

# Plot Kaplan-Meier curves
ggsurvplot(km_fit, conf.int = TRUE, pval = TRUE, risk.table = TRUE, censor.size = 0.6, size = 0.5,
           legend.labs = c("Depression", "Encephalopathy", "MCI", "Vascular"), legend.title = "Cluster",  
           palette = c("#26547c", "#ffd166", "#06d6a0", "#ef476f"), 
           title = "Kaplan-Meier Curve for AD to death Survival", 
           risk.table.height = .3)

# COMMAND ----------

AD_to_death = pat_cluster_unique_info %>% filter(cluster %in% c("Depression", "Encephalopathy")) %>%
  mutate(time = case_when(
    if_death == 1 ~ AD_to_death,
    if_death == 0 ~ AD_to_last_visit
  )) %>% mutate(time = as.numeric(time)) %>%
  select(PatientID, time, cluster, if_death) %>% 
  mutate(status = case_when(
    if_death == 1 ~ 1,
    if_death == 0 ~ 0
  )) 
# Fit Kaplan-Meier curves by groups
km_fit = survfit(Surv(time, status) ~ cluster, data = AD_to_death)
# Plot Kaplan-Meier curves
ggsurvplot(km_fit, conf.int = TRUE, pval = TRUE, risk.table = TRUE, censor.size = 0.6, size = 0.5,
           legend.labs = c("Depression", "Encephalopathy"), legend.title = "Cluster",  
           palette = c("dodgerblue2", "orchid2"), 
           title = "Kaplan-Meier Curve for AD to death Survival", 
           risk.table.height = .2)

# COMMAND ----------

AD_to_death = pat_cluster_unique_info %>% filter(cluster %in% c("Depression", "MCI")) %>%
  mutate(time = case_when(
    if_death == 1 ~ AD_to_death,
    if_death == 0 ~ AD_to_last_visit
  )) %>% mutate(time = as.numeric(time)) %>%
  select(PatientID, time, cluster, if_death) %>% 
  mutate(status = case_when(
    if_death == 1 ~ 1,
    if_death == 0 ~ 0
  )) 
# Fit Kaplan-Meier curves by groups
km_fit = survfit(Surv(time, status) ~ cluster, data = AD_to_death)
# Plot Kaplan-Meier curves
ggsurvplot(km_fit, conf.int = TRUE, pval = TRUE, risk.table = TRUE, censor.size = 0.6, size = 0.5,
           legend.labs = c("Depression", "MCI"), legend.title = "Cluster",  
           palette = c("dodgerblue2", "lightgreen"), 
           title = "Kaplan-Meier Curve for AD to death Survival", 
           risk.table.height = .2)

# COMMAND ----------

AD_to_death = pat_cluster_unique_info %>% filter(cluster %in% c("Depression", "Vascular")) %>%
  mutate(time = case_when(
    if_death == 1 ~ AD_to_death,
    if_death == 0 ~ AD_to_last_visit
  )) %>% mutate(time = as.numeric(time)) %>%
  select(PatientID, time, cluster, if_death) %>% 
  mutate(status = case_when(
    if_death == 1 ~ 1,
    if_death == 0 ~ 0
  )) 
# Fit Kaplan-Meier curves by groups
km_fit = survfit(Surv(time, status) ~ cluster, data = AD_to_death)
# Plot Kaplan-Meier curves
ggsurvplot(km_fit, conf.int = TRUE, pval = TRUE, risk.table = TRUE, censor.size = 0.6, size = 0.5,
           legend.labs = c("Depression", "Vascular"), legend.title = "Cluster",  
           palette = c("dodgerblue2", "orange"), 
           title = "Kaplan-Meier Curve for AD to death Survival", 
           risk.table.height = .2)

# COMMAND ----------

AD_to_death = pat_cluster_unique_info %>% filter(cluster %in% c("Encephalopathy", "MCI")) %>%
  mutate(time = case_when(
    if_death == 1 ~ AD_to_death,
    if_death == 0 ~ AD_to_last_visit
  )) %>% mutate(time = as.numeric(time)) %>%
  select(PatientID, time, cluster, if_death) %>% 
  mutate(status = case_when(
    if_death == 1 ~ 1,
    if_death == 0 ~ 0
  )) 
# Fit Kaplan-Meier curves by groups
km_fit = survfit(Surv(time, status) ~ cluster, data = AD_to_death)
# Plot Kaplan-Meier curves
ggsurvplot(km_fit, conf.int = TRUE, pval = TRUE, risk.table = TRUE, censor.size = 0.6, size = 0.5,
           legend.labs = c("Encephalopathy", "MCI"), legend.title = "Cluster",  
           palette = c("orchid2", "lightgreen"), 
           title = "Kaplan-Meier Curve for AD to death Survival", 
           risk.table.height = .2)

# COMMAND ----------

AD_to_death = pat_cluster_unique_info %>% filter(cluster %in% c("Encephalopathy", "Vascular")) %>%
  mutate(time = case_when(
    if_death == 1 ~ AD_to_death,
    if_death == 0 ~ AD_to_last_visit
  )) %>% mutate(time = as.numeric(time)) %>%
  select(PatientID, time, cluster, if_death) %>% 
  mutate(status = case_when(
    if_death == 1 ~ 1,
    if_death == 0 ~ 0
  )) 
# Fit Kaplan-Meier curves by groups
km_fit = survfit(Surv(time, status) ~ cluster, data = AD_to_death)
# Plot Kaplan-Meier curves
ggsurvplot(km_fit, conf.int = TRUE, pval = TRUE, risk.table = TRUE, censor.size = 0.6, size = 0.5,
           legend.labs = c("Encephalopathy", "Vascular"), legend.title = "Cluster",  
           palette = c("orchid2", "orange"), 
           title = "Kaplan-Meier Curve for AD to death Survival", 
           risk.table.height = .2)

# COMMAND ----------

AD_to_death = pat_cluster_unique_info %>% filter(cluster %in% c("MCI", "Vascular")) %>%
  mutate(time = case_when(
    if_death == 1 ~ AD_to_death,
    if_death == 0 ~ AD_to_last_visit
  )) %>% mutate(time = as.numeric(time)) %>%
  select(PatientID, time, cluster, if_death) %>% 
  mutate(status = case_when(
    if_death == 1 ~ 1,
    if_death == 0 ~ 0
  )) 
# Fit Kaplan-Meier curves by groups
km_fit = survfit(Surv(time, status) ~ cluster, data = AD_to_death)
# Plot Kaplan-Meier curves
ggsurvplot(km_fit, conf.int = TRUE, pval = TRUE, risk.table = TRUE, censor.size = 0.6, size = 0.5,
           legend.labs = c("MCI", "Vascular"), legend.title = "Cluster",  
           palette = c("lightgreen", "orange"), 
           title = "Kaplan-Meier Curve for AD to death Survival", 
           risk.table.height = .2)

# COMMAND ----------

AD_to_death = pat_cluster_unique_info %>% 
  mutate(time = case_when(
    if_death == 1 ~ AD_to_death,
    if_death == 0 ~ AD_to_last_visit
  )) %>% mutate(time = as.numeric(time)) %>%
  select(PatientID, time, cluster, if_death) %>% 
  mutate(status = case_when(
    if_death == 1 ~ 1,
    if_death == 0 ~ 0
  )) 
table(AD_to_death$status, AD_to_death$cluster)

# COMMAND ----------

# MAGIC %md
# MAGIC ### Diagnoses in deceased AD

# COMMAND ----------

# Function to perform hypothesis tests for categorical variables
categorical_hypothesis_tests = function(data, vars) {
  results = lapply(vars, function(var) {
    chisq_test = chisq.test(table(data[[var]], data[["if_death"]]))
    data.frame(
      variable = var,
      p_value = chisq_test$p.value
    )
  })
  do.call(rbind, results)
}

# COMMAND ----------

cluster_i = "Depression"

# COMMAND ----------

# for patients who have died, what diagnoses they have (by clusters)
cluster_death_pat = AD_to_death %>% filter(if_death == 1 & cluster == cluster_i) %>% pull(PatientID) %>% unique()
cluster_alive_pat = AD_to_death %>% filter(if_death == 0 & cluster == cluster_i) %>% pull(PatientID) %>% unique()
cluster_death_encounter = AD_encounters_filter %>% mutate(ICD_3digit = substr(ICD, 1, 3)) %>%
  select(-c(ICD, date_first_visit)) %>% unique() %>% group_by(PatientID, ICD_3digit) %>%
  arrange(EncounterDate) %>% slice(1) %>% ungroup() %>%
  filter(PatientID %in% c(cluster_death_pat, cluster_alive_pat)) %>% left_join(AD_demographic) %>% 
  mutate(if_death = case_when(
    PatientID %in% cluster_death_pat ~ 1,
    TRUE ~ 0
  )) %>% filter(EncounterDate <= AD_Date) %>% 
  select(PatientID, ICD_3digit, if_death) %>% unique() %>% mutate(status = 1) %>%
  pivot_wider(names_from = ICD_3digit, values_from = status, values_fill = 0) %>% 
  mutate_all(~ifelse(is.na(.), 0, .))
print(dim(cluster_death_encounter)) # dim = (880, 970)
# Calculate the sum of each column
column_sums = colSums(cluster_death_encounter[,-(1:2)])
# Print the sum of each column (only keep icds with prevalence >=10%)
select_icds = column_sums[column_sums >= 0.1*nrow(cluster_death_encounter)] %>% sort() %>% names() 
select_icds = setdiff(select_icds, "G30")
categorical_tests = categorical_hypothesis_tests(cluster_death_encounter, select_icds) %>%
  mutate(p_value = as.numeric(p_value)) %>% mutate(fdr_adj_p = p.adjust(p_value, method = "fdr")) 
for (i in 1:length(select_icds)) {
  icd = select_icds[i]
  frequency = table(cluster_death_encounter$if_death, cluster_death_encounter[[icd]])
  formatted_frequency = format(frequency, big.mark = ",")
  percentage = prop.table(frequency, margin = 1) * 100
  result_table = cbind(formatted_frequency[,2], percentage[,2]) %>% as.data.frame()
  names(result_table) = c("n", "percent")
  final_table = result_table %>% mutate(test_icd = paste0(n, " (", round(as.numeric(percent), 1), "%)")) %>%
    select(test_icd) %>% t() %>% as.data.frame() %>% rownames_to_column(var = "variable") %>%
    mutate(p_value = categorical_tests$fdr_adj_p[i]) %>%
    mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value,3)))) 
  names(final_table) = c("variable", "alive", "deceased", "adj_p_value")
  final_table$variable = icd
  if (i == 1) {
    death_summary_tbl = final_table
  } else {
    death_summary_tbl = rbind(death_summary_tbl, final_table)
  }
}
death_summary_tbl %>% filter(adj_p_value < 0.1) %>% arrange(variable)

# COMMAND ----------

cluster_i = "Encephalopathy"

# COMMAND ----------

# for patients who have died, what diagnoses they have (by clusters)
cluster_death_pat = AD_to_death %>% filter(if_death == 1 & cluster == cluster_i) %>% pull(PatientID) %>% unique()
cluster_alive_pat = AD_to_death %>% filter(if_death == 0 & cluster == cluster_i) %>% pull(PatientID) %>% unique()
cluster_death_encounter = AD_encounters_filter %>% mutate(ICD_3digit = substr(ICD, 1, 3)) %>%
  select(-c(ICD, date_first_visit)) %>% unique() %>% group_by(PatientID, ICD_3digit) %>%
  arrange(EncounterDate) %>% slice(1) %>% ungroup() %>%
  filter(PatientID %in% c(cluster_death_pat, cluster_alive_pat)) %>% left_join(AD_demographic) %>% 
  mutate(if_death = case_when(
    PatientID %in% cluster_death_pat ~ 1,
    TRUE ~ 0
  )) %>% filter(EncounterDate <= AD_Date) %>% 
  select(PatientID, ICD_3digit, if_death) %>% unique() %>% mutate(status = 1) %>%
  pivot_wider(names_from = ICD_3digit, values_from = status, values_fill = 0) %>% 
  mutate_all(~ifelse(is.na(.), 0, .))
print(dim(cluster_death_encounter)) # dim = (880, 970)
# Calculate the sum of each column
column_sums = colSums(cluster_death_encounter[,-(1:2)])
# Print the sum of each column (only keep icds with prevalence >=10%)
select_icds = column_sums[column_sums >= 0.1*nrow(cluster_death_encounter)] %>% sort() %>% names() 
select_icds = setdiff(select_icds, "G30")
categorical_tests = categorical_hypothesis_tests(cluster_death_encounter, select_icds) %>%
  mutate(p_value = as.numeric(p_value)) %>% mutate(fdr_adj_p = p.adjust(p_value, method = "fdr")) 
for (i in 1:length(select_icds)) {
  icd = select_icds[i]
  frequency = table(cluster_death_encounter$if_death, cluster_death_encounter[[icd]])
  formatted_frequency = format(frequency, big.mark = ",")
  percentage = prop.table(frequency, margin = 1) * 100
  result_table = cbind(formatted_frequency[,2], percentage[,2]) %>% as.data.frame()
  names(result_table) = c("n", "percent")
  final_table = result_table %>% mutate(test_icd = paste0(n, " (", round(as.numeric(percent), 1), "%)")) %>%
    select(test_icd) %>% t() %>% as.data.frame() %>% rownames_to_column(var = "variable") %>%
    mutate(p_value = categorical_tests$fdr_adj_p[i]) %>%
    mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value,3)))) 
  names(final_table) = c("variable", "alive", "deceased", "adj_p_value")
  final_table$variable = icd
  if (i == 1) {
    death_summary_tbl = final_table
  } else {
    death_summary_tbl = rbind(death_summary_tbl, final_table)
  }
}
death_summary_tbl %>% filter(adj_p_value < 0.1) %>% arrange(variable)

# COMMAND ----------

cluster_i = "MCI"

# COMMAND ----------

# for patients who have died, what diagnoses they have (by clusters)
cluster_death_pat = AD_to_death %>% filter(if_death == 1 & cluster == cluster_i) %>% pull(PatientID) %>% unique()
cluster_alive_pat = AD_to_death %>% filter(if_death == 0 & cluster == cluster_i) %>% pull(PatientID) %>% unique()
cluster_death_encounter = AD_encounters_filter %>% mutate(ICD_3digit = substr(ICD, 1, 3)) %>%
  select(-c(ICD, date_first_visit)) %>% unique() %>% group_by(PatientID, ICD_3digit) %>%
  arrange(EncounterDate) %>% slice(1) %>% ungroup() %>%
  filter(PatientID %in% c(cluster_death_pat, cluster_alive_pat)) %>% left_join(AD_demographic) %>% 
  mutate(if_death = case_when(
    PatientID %in% cluster_death_pat ~ 1,
    TRUE ~ 0
  )) %>% filter(EncounterDate <= AD_Date) %>% 
  select(PatientID, ICD_3digit, if_death) %>% unique() %>% mutate(status = 1) %>%
  pivot_wider(names_from = ICD_3digit, values_from = status, values_fill = 0) %>% 
  mutate_all(~ifelse(is.na(.), 0, .))
print(dim(cluster_death_encounter)) # dim = (880, 970)
# Calculate the sum of each column
column_sums = colSums(cluster_death_encounter[,-(1:2)])
# Print the sum of each column (only keep icds with prevalence >=10%)
select_icds = column_sums[column_sums >= 0.1*nrow(cluster_death_encounter)] %>% sort() %>% names() 
select_icds = setdiff(select_icds, "G30")
categorical_tests = categorical_hypothesis_tests(cluster_death_encounter, select_icds) %>%
  mutate(p_value = as.numeric(p_value)) %>% mutate(fdr_adj_p = p.adjust(p_value, method = "fdr")) 
for (i in 1:length(select_icds)) {
  icd = select_icds[i]
  frequency = table(cluster_death_encounter$if_death, cluster_death_encounter[[icd]])
  formatted_frequency = format(frequency, big.mark = ",")
  percentage = prop.table(frequency, margin = 1) * 100
  result_table = cbind(formatted_frequency[,2], percentage[,2]) %>% as.data.frame()
  names(result_table) = c("n", "percent")
  final_table = result_table %>% mutate(test_icd = paste0(n, " (", round(as.numeric(percent), 1), "%)")) %>%
    select(test_icd) %>% t() %>% as.data.frame() %>% rownames_to_column(var = "variable") %>%
    mutate(p_value = categorical_tests$fdr_adj_p[i]) %>%
    mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value,3)))) 
  names(final_table) = c("variable", "alive", "deceased", "adj_p_value")
  final_table$variable = icd
  if (i == 1) {
    death_summary_tbl = final_table
  } else {
    death_summary_tbl = rbind(death_summary_tbl, final_table)
  }
}
death_summary_tbl %>% filter(adj_p_value < 0.1) %>% arrange(variable)

# COMMAND ----------

cluster_i = "Vascular"

# COMMAND ----------

# for patients who have died, what diagnoses they have (by clusters)
cluster_death_pat = AD_to_death %>% filter(if_death == 1 & cluster == cluster_i) %>% pull(PatientID) %>% unique()
cluster_alive_pat = AD_to_death %>% filter(if_death == 0 & cluster == cluster_i) %>% pull(PatientID) %>% unique()
cluster_death_encounter = AD_encounters_filter %>% mutate(ICD_3digit = substr(ICD, 1, 3)) %>%
  select(-c(ICD, date_first_visit)) %>% unique() %>% group_by(PatientID, ICD_3digit) %>%
  arrange(EncounterDate) %>% slice(1) %>% ungroup() %>%
  filter(PatientID %in% c(cluster_death_pat, cluster_alive_pat)) %>% left_join(AD_demographic) %>% 
  mutate(if_death = case_when(
    PatientID %in% cluster_death_pat ~ 1,
    TRUE ~ 0
  )) %>% filter(EncounterDate <= AD_Date) %>% 
  select(PatientID, ICD_3digit, if_death) %>% unique() %>% mutate(status = 1) %>%
  pivot_wider(names_from = ICD_3digit, values_from = status, values_fill = 0) %>% 
  mutate_all(~ifelse(is.na(.), 0, .))
print(dim(cluster_death_encounter)) # dim = (880, 970)
# Calculate the sum of each column
column_sums = colSums(cluster_death_encounter[,-(1:2)])
# Print the sum of each column (only keep icds with prevalence >=10%)
select_icds = column_sums[column_sums >= 0.1*nrow(cluster_death_encounter)] %>% sort() %>% names() 
select_icds = setdiff(select_icds, "G30")
categorical_tests = categorical_hypothesis_tests(cluster_death_encounter, select_icds) %>%
  mutate(p_value = as.numeric(p_value)) %>% mutate(fdr_adj_p = p.adjust(p_value, method = "fdr")) 
for (i in 1:length(select_icds)) {
  icd = select_icds[i]
  frequency = table(cluster_death_encounter$if_death, cluster_death_encounter[[icd]])
  formatted_frequency = format(frequency, big.mark = ",")
  percentage = prop.table(frequency, margin = 1) * 100
  result_table = cbind(formatted_frequency[,2], percentage[,2]) %>% as.data.frame()
  names(result_table) = c("n", "percent")
  final_table = result_table %>% mutate(test_icd = paste0(n, " (", round(as.numeric(percent), 1), "%)")) %>%
    select(test_icd) %>% t() %>% as.data.frame() %>% rownames_to_column(var = "variable") %>%
    mutate(p_value = categorical_tests$fdr_adj_p[i]) %>%
    mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value,3)))) 
  names(final_table) = c("variable", "alive", "deceased", "adj_p_value")
  final_table$variable = icd
  if (i == 1) {
    death_summary_tbl = final_table
  } else {
    death_summary_tbl = rbind(death_summary_tbl, final_table)
  }
}
death_summary_tbl %>% filter(adj_p_value < 0.1) %>% arrange(variable)

# COMMAND ----------

# MAGIC %md
# MAGIC ### Box plots from central node to AD

# COMMAND ----------

library(remotes)
# install_version("Matrix", version = "1.6-5")
install.packages("ggpubr")
library(ggpubr)

# COMMAND ----------

AD_encounters_3digit = AD_encounters_filter %>% mutate(ICD_3digit = substr(ICD, 1, 3)) %>%
  select(-c(ICD, date_first_visit)) %>% unique() %>% group_by(PatientID, ICD_3digit) %>%
  arrange(EncounterDate) %>% slice(1) %>% ungroup() %>%
  select(PatientID, EncounterDate, ICD_3digit, AD_Date) %>% unique()
print(dim(AD_encounters_3digit)) # dim = (729347,4)

# COMMAND ----------

central_nodes = c("F32", "I67", "G93", "G31")

# COMMAND ----------

pat_cluster_unique = pat_traj_final_cluster %>% filter(PatientID %!in% overlap_pat)
for (i in 1:4) {
  cluster_i = i
  pat_cluster_id = pat_cluster_unique %>% filter(Kmeans == cluster_i) %>% 
    pull(PatientID) %>% unique() 
  # print(length(pat_cluster_id))
  central_encounter = AD_encounters_3digit %>% filter(PatientID %in% pat_cluster_id) %>%
    filter(ICD_3digit == central_nodes[i]) %>% 
    mutate(central_to_AD = as.numeric(AD_Date - EncounterDate)/365.25) %>%
    select(PatientID, central_to_AD) %>% mutate(central_node = central_nodes[i]) %>% unique()
  if (i == 1) {
    central_time_info = central_encounter
  } else {
    central_time_info = rbind(central_time_info, central_encounter)
  }
}
print(dim(central_time_info)) # dim = (1626,3)
table(central_time_info$central_node)

# COMMAND ----------

central_time_info %>% group_by(central_node) %>%
  summarise(mean_time = mean(central_to_AD, na.rm = T),
            median_time = median(central_to_AD, na.rm = T),
            max_time = max(central_to_AD, na.rm = T),
            min_time = min(central_to_AD, na.rm = T))

# COMMAND ----------

head(central_time_info)

# COMMAND ----------

central_time_info$central_node = factor(central_time_info$central_node, levels = c("F32", "G93", "G31", "I67")) # Specify the order of categories

# COMMAND ----------

# Perorm pairwise comparisons
compare_means(central_to_AD ~ central_node,  data = central_time_info)

# COMMAND ----------

# Visualize: Specify the comparisons you want
my_comparisons = list(c("F32", "I67"), c("F32", "G93"), c("F32", "G31"),
                      c("I67", "G93"), c("I67", "G31"), c("G93", "G31") )
central_time_info %>%
  ggplot(aes(x = central_node, y = central_to_AD, fill = central_node)) + 
    geom_boxplot(alpha = 0.5) + 
    scale_fill_manual(values = c("F32" = "#26547c", "G93" = "#ffd166", "G31" = "#06d6a0", "I67" = "#ef476f")) +
    theme(legend.position = "none", panel.background = element_blank(),  # Remove panel background
    plot.background = element_blank(),   # Remove plot background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(),          # Add axis lines
    axis.ticks = element_line(),         # Add axis ticks
    axis.title = element_text(),         # Customize axis title
    axis.text = element_text()) + 
    xlab("") + ylab("Time from central node to AD") +
    stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
    stat_compare_means(label.y = 18)     # Add global p-value

# COMMAND ----------

# MAGIC %md
# MAGIC # 4. Multi-cluster patients

# COMMAND ----------

# MAGIC %md
# MAGIC ### MCI cluster

# COMMAND ----------

mci_cluster = pat_traj_final_cluster %>% filter(PatientID %in% pat_cluster_4) %>% 
  mutate(cluster = case_when(
    Kmeans == 1 ~ "Depression",
    Kmeans == 2 ~ "Vascular",
    Kmeans == 3 ~ "Encephalopathy",
    Kmeans == 4 ~ "MCI"
  )) %>% 
  select(PatientID, cluster) %>% mutate(value = 1) %>% unique() %>%
  pivot_wider(names_from = cluster,
    values_from = value,
    values_fill = list(value = 0)) %>%
  mutate(group = case_when(
    Vascular == 0 & Depression == 0 & Encephalopathy == 0 ~ "MCI_only",
    Vascular == 1 & Depression == 0 & Encephalopathy == 0 ~ "MCI_vascular",
    Vascular == 0 & Depression == 1 & Encephalopathy == 0 ~ "MCI_depression",
    Vascular == 0 & Depression == 0 & Encephalopathy == 1 ~ "MCI_encephalopathy",
    Vascular == 1 & Depression == 1 & Encephalopathy == 0 ~ "MCI_vascular_depression",
    Vascular == 1 & Depression == 0 & Encephalopathy == 1 ~ "MCI_vascular_encephalopathy",
    Vascular == 0 & Depression == 1 & Encephalopathy == 1 ~ "MCI_depression_encephalopathy",
  )) %>% select(PatientID, group) %>% unique() %>%
  mutate(group_general = if_else(group == "MCI_only", "MCI_only", "MCI_multi"))
print(dim(mci_cluster)) # dim = (1502,2)
head(mci_cluster)

# COMMAND ----------

table(mci_cluster$group)

# COMMAND ----------

table(mci_cluster$group_general)

# COMMAND ----------

mci_pat_demo = mci_cluster %>% left_join(AD_demographic) %>%
  mutate(age_AD = as.numeric(AD_Date - BirthDate)/365.25) %>%
  mutate(AD_to_last = as.numeric(date_last_visit - AD_Date)/365.25) %>%
  mutate(AD_to_death = as.numeric(DeathDate - AD_Date)/365.25) %>%
  mutate(first_visit_to_AD = as.numeric(AD_Date - date_first_visit)/365.25) %>%
  select(PatientID, female, race_ethnicity, PatientLivingStatus, n_icd_3digit, n_encounter, record_length, 
        enc_per_yr, age_last_visit, location_source_value, first_visit_to_AD, age_AD, AD_to_last, AD_to_death, 
        group, group_general)
print(dim(mci_pat_demo)) # dim = (1502,16)

# COMMAND ----------

# List of columns to compare
continuous_vars = c("age_last_visit", "record_length", "n_encounter", "n_icd_3digit", "enc_per_yr")
# Function to perform hypothesis tests for continuous variables
continuous_hypothesis_tests = function(data, vars) {
  results = lapply(vars, function(var) {
    wilcox_test = kruskal.test(data[[var]] ~ data[["group"]])
    data.frame(
      variable = var,
      p_value = wilcox_test$p.value
    )
  })
  do.call(rbind, results)
}
categorical_vars = c("female", "race_ethnicity", "PatientLivingStatus", "location_source_value")
# Function to summarize categorical variables
# Function to perform hypothesis tests for categorical variables
categorical_hypothesis_tests = function(data, vars) {
  results = lapply(vars, function(var) {
    chisq_test = chisq.test(table(data[[var]], data[["group"]]))
    data.frame(
      variable = var,
      p_value = chisq_test$p.value
    )
  })
  do.call(rbind, results)
}

# COMMAND ----------

#======== Continuous variables =============
cont_median = mci_pat_demo %>% group_by(group) %>%
    summarise(across(all_of(continuous_vars), list(median = median)), .groups = "drop") %>% t() %>%as.data.frame() 
colnames(cont_median) = cont_median[1, ]
cont_median = cont_median[-1, ]
row.names(cont_median) = NULL
cont_median = cont_median %>% 
  mutate(MCI_depression_m = round(as.numeric(MCI_depression), 1), MCI_depression_encephalopathy_m = round(as.numeric(MCI_depression_encephalopathy), 1), 
        MCI_encephalopathy_m = round(as.numeric(MCI_encephalopathy), 1), MCI_only_m = round(as.numeric(MCI_only), 1), 
        MCI_vascular_m = round(as.numeric(MCI_vascular), 1), MCI_vascular_depression_m = round(as.numeric(MCI_vascular_depression), 1), 
        MCI_vascular_encephalopathy_m = round(as.numeric(MCI_vascular_encephalopathy), 1), variable = continuous_vars) %>%
  select(variable, MCI_only_m, MCI_depression_m, MCI_encephalopathy_m, MCI_vascular_m, 
        MCI_depression_encephalopathy_m, MCI_vascular_depression_m, MCI_vascular_encephalopathy_m)

cont_IQR = mci_pat_demo %>% group_by(group) %>%
    summarise(across(all_of(continuous_vars), list(IQR = IQR)), .groups = "drop") %>% t() %>%as.data.frame() 
colnames(cont_IQR) = cont_IQR[1, ]
cont_IQR = cont_IQR[-1, ]
row.names(cont_IQR) = NULL
cont_IQR = cont_IQR %>% 
  mutate(MCI_depression_i = round(as.numeric(MCI_depression), 1), MCI_depression_encephalopathy_i = round(as.numeric(MCI_depression_encephalopathy), 1), 
        MCI_encephalopathy_i = round(as.numeric(MCI_encephalopathy), 1), MCI_only_i = round(as.numeric(MCI_only), 1), 
        MCI_vascular_i = round(as.numeric(MCI_vascular), 1), MCI_vascular_depression_i = round(as.numeric(MCI_vascular_depression), 1), 
        MCI_vascular_encephalopathy_i = round(as.numeric(MCI_vascular_encephalopathy), 1), variable = continuous_vars) %>%
  select(variable, MCI_only_i, MCI_depression_i, MCI_encephalopathy_i, MCI_vascular_i, 
        MCI_depression_encephalopathy_i, MCI_vascular_depression_i, MCI_vascular_encephalopathy_i)

cont_tests = continuous_hypothesis_tests(mci_pat_demo, continuous_vars)

cont_full = cont_median %>% full_join(cont_IQR) %>%
  mutate(MCI_depression = paste0(MCI_depression_m, " (", MCI_depression_i, ")"), 
        MCI_depression_encephalopathy = paste0(MCI_depression_encephalopathy_m, " (", MCI_depression_encephalopathy_i, ")"),
        MCI_encephalopathy = paste0(MCI_encephalopathy_m, " (", MCI_encephalopathy_i, ")"),
        MCI_only = paste0(MCI_only_m, " (", MCI_only_i, ")"),
        MCI_vascular = paste0(MCI_vascular_m, " (", MCI_vascular_i, ")"),
        MCI_vascular_depression = paste0(MCI_vascular_depression_m, " (", MCI_vascular_depression_i, ")"), 
        MCI_vascular_encephalopathy = paste0(MCI_vascular_encephalopathy_m, " (", MCI_vascular_encephalopathy_i, ")")) %>%
  select(variable, MCI_only, MCI_depression, MCI_encephalopathy, MCI_vascular, 
        MCI_depression_encephalopathy, MCI_vascular_depression, MCI_vascular_encephalopathy) %>% full_join(cont_tests) %>%
  mutate(p_value = as.numeric(p_value)) %>%
  mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value, 3)))) 

# COMMAND ----------

#======== Categorical variables =============
categorical_tests = categorical_hypothesis_tests(mci_pat_demo, categorical_vars)
categorical_tests = categorical_tests %>%
  mutate(p_value = as.numeric(p_value)) %>%
  mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value,3)))) 

# Gender
frequency = table(mci_pat_demo$group, mci_pat_demo$female)
formatted_frequency = format(frequency, big.mark = ",")
percentage = prop.table(frequency, margin = 1) * 100
result_table = cbind(formatted_frequency[,2], percentage[,2]) %>% as.data.frame()
names(result_table) = c("n", "percent")
gender_table = result_table %>% mutate(female = paste0(n, " (", round(as.numeric(percent), 1), "%)")) %>%
  select(female) %>% t() %>% as.data.frame() %>% rownames_to_column(var = "variable") %>%
  mutate(p_value = categorical_tests$p_value[1])
names(gender_table) = c("variable", "MCI_depression", "MCI_depression_encephalopathy", "MCI_encephalopathy", "MCI_only", "MCI_vascular", 
        "MCI_vascular_depression", "MCI_vascular_encephalopathy", "p_value")

# Race-ethnicity
frequency = table(mci_pat_demo$group, mci_pat_demo$race_ethnicity) 
formatted_frequency = format(frequency, big.mark = ",") %>% as.data.frame()
percentage = prop.table(frequency, margin = 1) * 100
formatted_percentage = format(round(percentage, 1)) %>% as.data.frame()
result_table = data.frame()
col_names = names(formatted_percentage)
for (i in 1:(dim(frequency)[2])) {
  if (i == 1) {result_table = cbind(formatted_frequency[,i], formatted_percentage[,i]) %>% as.data.frame()}
  else (result_table = cbind(result_table, formatted_frequency[,i], formatted_percentage[,i]) %>% as.data.frame())
  result_table[[col_names[i]]] = paste0(trimws(result_table[,ncol(result_table)-1]), " (", trimws(result_table[,ncol(result_table)]), "%)")
}
# Select every third column starting from the third column
race_table = result_table[, seq(3, ncol(result_table), by = 3)] %>% t() %>% as.data.frame()
names(race_table) = c("MCI_depression", "MCI_depression_encephalopathy", "MCI_encephalopathy", "MCI_only", "MCI_vascular", 
        "MCI_vascular_depression", "MCI_vascular_encephalopathy")
race_table = race_table %>% mutate(p_value = categorical_tests$p_value[2]) %>% rownames_to_column(var = "variable")

# Death
frequency = table(mci_pat_demo$group, mci_pat_demo$PatientLivingStatus)
formatted_frequency = format(frequency, big.mark = ",")
percentage = prop.table(frequency, margin = 1) * 100
result_table = cbind(formatted_frequency[,2], percentage[,2]) %>% as.data.frame()
names(result_table) = c("n", "percent")
death_table = result_table %>% mutate(death = paste0(n, " (", round(as.numeric(percent), 1), "%)")) %>%
  select(death) %>% t() %>% as.data.frame() %>% rownames_to_column(var = "variable") %>%
  mutate(p_value = categorical_tests$p_value[3])
names(death_table) = c("variable", "MCI_depression", "MCI_depression_encephalopathy", "MCI_encephalopathy", "MCI_only", "MCI_vascular", 
        "MCI_vascular_depression", "MCI_vascular_encephalopathy", "p_value")
# location
frequency = table(mci_pat_demo$group, mci_pat_demo$location_source_value) 
formatted_frequency = format(frequency, big.mark = ",") %>% as.data.frame()
percentage = prop.table(frequency, margin = 1) * 100
formatted_percentage = format(round(percentage, 1)) %>% as.data.frame()
result_table = data.frame()
col_names = names(formatted_percentage)
for (i in 1:(dim(frequency)[2])) {
  if (i == 1) {result_table = cbind(formatted_frequency[,i], formatted_percentage[,i]) %>% as.data.frame()}
  else (result_table = cbind(result_table, formatted_frequency[,i], formatted_percentage[,i]) %>% as.data.frame())
  result_table[[col_names[i]]] = paste0(trimws(result_table[,ncol(result_table)-1]), " (", trimws(result_table[,ncol(result_table)]), "%)")
}
# Select every third column starting from the third column
location_table = result_table[, seq(3, ncol(result_table), by = 3)] %>% t() %>% as.data.frame()
names(location_table) = c("MCI_depression", "MCI_depression_encephalopathy", "MCI_encephalopathy", "MCI_only", "MCI_vascular", 
        "MCI_vascular_depression", "MCI_vascular_encephalopathy")
location_table = location_table %>% mutate(p_value = categorical_tests$p_value[4]) %>% rownames_to_column(var = "variable")

# COMMAND ----------

final_sumstat_unique = rbind(cont_full, gender_table, death_table, c("Race-ethncity", "", "", "", "", "", "", "", ""), 
                race_table, c("Location sources", "", "", "", "", "", "", "", ""), location_table)
final_sumstat_unique

# COMMAND ----------

# List of columns to compare
continuous_vars = c("age_AD", "first_visit_to_AD", "AD_to_last", "AD_to_death")
# Function to perform hypothesis tests for continuous variables
continuous_hypothesis_tests = function(data, vars) {
  results = lapply(vars, function(var) {
    wilcox_test = kruskal.test(data[[var]] ~ data[["group"]])
    data.frame(
      variable = var,
      p_value = wilcox_test$p.value
    )
  })
  do.call(rbind, results)
}
#======== Continuous variables =============
cont_median = mci_pat_demo %>% group_by(group) %>%
    summarise(across(all_of(continuous_vars), list(median = median), na.rm = T), .groups = "drop") %>% t() %>%as.data.frame() 
colnames(cont_median) = cont_median[1, ]
cont_median = cont_median[-1, ]
row.names(cont_median) = NULL
cont_median = cont_median %>% 
  mutate(MCI_depression_m = round(as.numeric(MCI_depression), 1), MCI_depression_encephalopathy_m = round(as.numeric(MCI_depression_encephalopathy), 1), 
        MCI_encephalopathy_m = round(as.numeric(MCI_encephalopathy), 1), MCI_only_m = round(as.numeric(MCI_only), 1), 
        MCI_vascular_m = round(as.numeric(MCI_vascular), 1), MCI_vascular_depression_m = round(as.numeric(MCI_vascular_depression), 1), 
        MCI_vascular_encephalopathy_m = round(as.numeric(MCI_vascular_encephalopathy), 1), variable = continuous_vars) %>%
  select(variable, MCI_only_m, MCI_depression_m, MCI_encephalopathy_m, MCI_vascular_m, 
        MCI_depression_encephalopathy_m, MCI_vascular_depression_m, MCI_vascular_encephalopathy_m)


cont_IQR = mci_pat_demo %>% group_by(group) %>%
    summarise(across(all_of(continuous_vars), list(IQR = IQR), na.rm = T), .groups = "drop") %>% t() %>%as.data.frame() 
colnames(cont_IQR) = cont_IQR[1, ]
cont_IQR = cont_IQR[-1, ]
row.names(cont_IQR) = NULL
cont_IQR = cont_IQR %>% 
  mutate(MCI_depression_i = round(as.numeric(MCI_depression), 1), MCI_depression_encephalopathy_i = round(as.numeric(MCI_depression_encephalopathy), 1), 
        MCI_encephalopathy_i = round(as.numeric(MCI_encephalopathy), 1), MCI_only_i = round(as.numeric(MCI_only), 1), 
        MCI_vascular_i = round(as.numeric(MCI_vascular), 1), MCI_vascular_depression_i = round(as.numeric(MCI_vascular_depression), 1), 
        MCI_vascular_encephalopathy_i = round(as.numeric(MCI_vascular_encephalopathy), 1), variable = continuous_vars) %>%
  select(variable, MCI_only_i, MCI_depression_i, MCI_encephalopathy_i, MCI_vascular_i, 
        MCI_depression_encephalopathy_i, MCI_vascular_depression_i, MCI_vascular_encephalopathy_i)

cont_tests = continuous_hypothesis_tests(mci_pat_demo, continuous_vars)
cont_full = cont_median %>% full_join(cont_IQR) %>%
  mutate(MCI_depression = paste0(MCI_depression_m, " (", MCI_depression_i, ")"), 
        MCI_depression_encephalopathy = paste0(MCI_depression_encephalopathy_m, " (", MCI_depression_encephalopathy_i, ")"),
        MCI_encephalopathy = paste0(MCI_encephalopathy_m, " (", MCI_encephalopathy_i, ")"),
        MCI_only = paste0(MCI_only_m, " (", MCI_only_i, ")"),
        MCI_vascular = paste0(MCI_vascular_m, " (", MCI_vascular_i, ")"),
        MCI_vascular_depression = paste0(MCI_vascular_depression_m, " (", MCI_vascular_depression_i, ")"), 
        MCI_vascular_encephalopathy = paste0(MCI_vascular_encephalopathy_m, " (", MCI_vascular_encephalopathy_i, ")")) %>%
  select(variable, MCI_only, MCI_depression, MCI_encephalopathy, MCI_vascular, 
        MCI_depression_encephalopathy, MCI_vascular_depression, MCI_vascular_encephalopathy) %>% full_join(cont_tests) %>%
  mutate(p_value = as.numeric(p_value)) %>%
  mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value, 3)))) 
cont_full

# COMMAND ----------

# List of columns to compare
continuous_vars = c("age_last_visit", "record_length", "n_encounter", "n_icd_3digit", "enc_per_yr")
# Function to perform hypothesis tests for continuous variables
continuous_hypothesis_tests = function(data, vars) {
  results = lapply(vars, function(var) {
    wilcox_test = kruskal.test(data[[var]] ~ data[["group_general"]])
    data.frame(
      variable = var,
      p_value = wilcox_test$p.value
    )
  })
  do.call(rbind, results)
}

# COMMAND ----------

#======== Continuous variables =============
cont_median = mci_pat_demo %>% group_by(group_general) %>%
    summarise(across(all_of(continuous_vars), list(median = median)), .groups = "drop") %>% t() %>%as.data.frame() 
colnames(cont_median) = cont_median[1, ]
cont_median = cont_median[-1, ]
row.names(cont_median) = NULL
cont_median = cont_median %>% 
  mutate(MCI_only_m = round(as.numeric(MCI_only), 1), MCI_multi_m = round(as.numeric(MCI_multi), 1), variable = continuous_vars) %>%
  select(variable, MCI_only_m, MCI_multi_m)

cont_IQR = mci_pat_demo %>% group_by(group_general) %>%
    summarise(across(all_of(continuous_vars), list(IQR = IQR)), .groups = "drop") %>% t() %>%as.data.frame() 
colnames(cont_IQR) = cont_IQR[1, ]
cont_IQR = cont_IQR[-1, ]
row.names(cont_IQR) = NULL
cont_IQR = cont_IQR %>% 
  mutate(MCI_only_i = round(as.numeric(MCI_only), 1), 
        MCI_multi_i = round(as.numeric(MCI_multi), 1), variable = continuous_vars) %>%
  select(variable, MCI_only_i, MCI_multi_i)

cont_tests = continuous_hypothesis_tests(mci_pat_demo, continuous_vars)

cont_full = cont_median %>% full_join(cont_IQR) %>%
  mutate(MCI_only = paste0(MCI_only_m, " (", MCI_only_i, ")"),
        MCI_multi = paste0(MCI_multi_m, " (", MCI_multi_i, ")")) %>%
  select(variable, MCI_only, MCI_multi) %>% full_join(cont_tests) %>%
  mutate(p_value = as.numeric(p_value)) %>%
  mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value, 3)))) 

# COMMAND ----------

cont_full

# COMMAND ----------

table(mci_pat_demo$group_general)

# COMMAND ----------

table(mci_pat_demo$group_general) %>% prop.table()

# COMMAND ----------

categorical_vars = c("female", "race_ethnicity", "PatientLivingStatus", "location_source_value")
# Function to summarize categorical variables
# Function to perform hypothesis tests for categorical variables
categorical_hypothesis_tests = function(data, vars) {
  results = lapply(vars, function(var) {
    chisq_test = chisq.test(table(data[[var]], data[["group_general"]]))
    data.frame(
      variable = var,
      p_value = chisq_test$p.value
    )
  })
  do.call(rbind, results)
}

# COMMAND ----------

#======== Categorical variables =============
categorical_tests = categorical_hypothesis_tests(mci_pat_demo, categorical_vars)
categorical_tests = categorical_tests %>%
  mutate(p_value = as.numeric(p_value)) %>%
  mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value,3)))) 

# Gender
frequency = table(mci_pat_demo$group_general, mci_pat_demo$female)
formatted_frequency = format(frequency, big.mark = ",")
percentage = prop.table(frequency, margin = 1) * 100
result_table = cbind(formatted_frequency[,2], percentage[,2]) %>% as.data.frame()
names(result_table) = c("n", "percent")
gender_table = result_table %>% mutate(female = paste0(n, " (", round(as.numeric(percent), 1), "%)")) %>%
  select(female) %>% t() %>% as.data.frame() %>% rownames_to_column(var = "variable") %>%
  mutate(p_value = categorical_tests$p_value[1])
names(gender_table) = c("variable", "MCI_multi", "MCI_only", "p_value")

# Race-ethnicity
frequency = table(mci_pat_demo$group_general, mci_pat_demo$race_ethnicity) 
formatted_frequency = format(frequency, big.mark = ",") %>% as.data.frame()
percentage = prop.table(frequency, margin = 1) * 100
formatted_percentage = format(round(percentage, 1)) %>% as.data.frame()
result_table = data.frame()
col_names = names(formatted_percentage)
for (i in 1:(dim(frequency)[2])) {
  if (i == 1) {result_table = cbind(formatted_frequency[,i], formatted_percentage[,i]) %>% as.data.frame()}
  else (result_table = cbind(result_table, formatted_frequency[,i], formatted_percentage[,i]) %>% as.data.frame())
  result_table[[col_names[i]]] = paste0(trimws(result_table[,ncol(result_table)-1]), " (", trimws(result_table[,ncol(result_table)]), "%)")
}
# Select every third column starting from the third column
race_table = result_table[, seq(3, ncol(result_table), by = 3)] %>% t() %>% as.data.frame()
names(race_table) = c("MCI_multi", "MCI_only")
race_table = race_table %>% mutate(p_value = categorical_tests$p_value[2]) %>% rownames_to_column(var = "variable")
# Death
frequency = table(mci_pat_demo$group_general, mci_pat_demo$PatientLivingStatus)
formatted_frequency = format(frequency, big.mark = ",")
percentage = prop.table(frequency, margin = 1) * 100
result_table = cbind(formatted_frequency[,2], percentage[,2]) %>% as.data.frame()
names(result_table) = c("n", "percent")
death_table = result_table %>% mutate(death = paste0(n, " (", round(as.numeric(percent), 1), "%)")) %>%
  select(death) %>% t() %>% as.data.frame() %>% rownames_to_column(var = "variable") %>%
  mutate(p_value = categorical_tests$p_value[3])
names(death_table) = c("variable", "MCI_multi", "MCI_only", "p_value")
# location
frequency = table(mci_pat_demo$group_general, mci_pat_demo$location_source_value) 
formatted_frequency = format(frequency, big.mark = ",") %>% as.data.frame()
percentage = prop.table(frequency, margin = 1) * 100
formatted_percentage = format(round(percentage, 1)) %>% as.data.frame()
result_table = data.frame()
col_names = names(formatted_percentage)
for (i in 1:(dim(frequency)[2])) {
  if (i == 1) {result_table = cbind(formatted_frequency[,i], formatted_percentage[,i]) %>% as.data.frame()}
  else (result_table = cbind(result_table, formatted_frequency[,i], formatted_percentage[,i]) %>% as.data.frame())
  result_table[[col_names[i]]] = paste0(trimws(result_table[,ncol(result_table)-1]), " (", trimws(result_table[,ncol(result_table)]), "%)")
}
# Select every third column starting from the third column
location_table = result_table[, seq(3, ncol(result_table), by = 3)] %>% t() %>% as.data.frame()
names(location_table) = c("MCI_multi", "MCI_only")
location_table = location_table %>% mutate(p_value = categorical_tests$p_value[4]) %>% rownames_to_column(var = "variable")

# COMMAND ----------

final_sumstat_unique = rbind(cont_full, gender_table, death_table, c("Race-ethncity", "", "", "", "", "", "", "", ""), 
                race_table, c("Location sources", "", "", "", "", "", "", "", ""), location_table)
final_sumstat_unique

# COMMAND ----------

# List of columns to compare
continuous_vars = c("age_AD", "first_visit_to_AD", "AD_to_last", "AD_to_death")
# Function to perform hypothesis tests for continuous variables
continuous_hypothesis_tests = function(data, vars) {
  results = lapply(vars, function(var) {
    wilcox_test = kruskal.test(data[[var]] ~ data[["group_general"]])
    data.frame(
      variable = var,
      p_value = wilcox_test$p.value
    )
  })
  do.call(rbind, results)
}
#======== Continuous variables =============
cont_median = mci_pat_demo %>% group_by(group_general) %>%
    summarise(across(all_of(continuous_vars), list(median = median), na.rm = T), .groups = "drop") %>% t() %>%as.data.frame() 
colnames(cont_median) = cont_median[1, ]
cont_median = cont_median[-1, ]
row.names(cont_median) = NULL
cont_median = cont_median %>% 
  mutate(MCI_only_m = round(as.numeric(MCI_only), 1), 
        MCI_multi_m = round(as.numeric(MCI_multi), 1), variable = continuous_vars) %>%
  select(variable, MCI_only_m, MCI_multi_m)


cont_IQR = mci_pat_demo %>% group_by(group_general) %>%
    summarise(across(all_of(continuous_vars), list(IQR = IQR), na.rm = T), .groups = "drop") %>% t() %>%as.data.frame() 
colnames(cont_IQR) = cont_IQR[1, ]
cont_IQR = cont_IQR[-1, ]
row.names(cont_IQR) = NULL
cont_IQR = cont_IQR %>% 
  mutate(MCI_only_i = round(as.numeric(MCI_only), 1), 
        MCI_multi_i = round(as.numeric(MCI_multi), 1), variable = continuous_vars) %>%
  select(variable, MCI_only_i, MCI_multi_i)

cont_tests = continuous_hypothesis_tests(mci_pat_demo, continuous_vars)
cont_full = cont_median %>% full_join(cont_IQR) %>%
  mutate(MCI_only = paste0(MCI_only_m, " (", MCI_only_i, ")"),
        MCI_multi = paste0(MCI_multi_m, " (", MCI_multi_i, ")")) %>%
  select(variable, MCI_only, MCI_multi) %>% full_join(cont_tests) %>%
  mutate(p_value = as.numeric(p_value)) %>%
  mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value, 3)))) 
cont_full

# COMMAND ----------

