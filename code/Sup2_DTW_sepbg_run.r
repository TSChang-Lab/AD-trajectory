# Databricks notebook source
# MAGIC %md
# MAGIC # Supplementary Code 2. DTW alignment of patient trajectories

# COMMAND ----------

# calculate the divided point
n = 3272
comparisons = (n * (n - 1)) / 2
comparisons / 5
# 350|400|470|610|1442

# COMMAND ----------

# MAGIC %md
# MAGIC ## Run in sep files

# COMMAND ----------

library(tidyverse)
raw_data_path = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/AD_trajectory/data/"
load(file = paste0(raw_data_path, "patient_data/final/patient_traj_final_0313.rda"))

# COMMAND ----------

load(file = paste0(raw_data_path, "upload_ucla/embed_sim_cal.rda"))
print(dim(embed_sim_cal)) # dim = (1619,1619)
embed_sim_cal[upper.tri(embed_sim_cal)] = t(embed_sim_cal)[upper.tri(embed_sim_cal)]
embed.dist = 1 - embed_sim_cal
embed_icd = rownames(embed_sim_cal) # N = 1619
# prepare unique trajectories
patient_traj_unique = patient_traj_final %>% select(-PatientID) %>% unique() 
print(dim(patient_traj_unique)) # dim = (3272,11)
# sort patient IDs
dtw_traj_long_rmdup = patient_traj_unique %>% 
  pivot_longer(cols = setdiff(names(patient_traj_unique), "traj_id"), names_to = "step", values_to = "ICD_3digit") %>% drop_na() %>% mutate(step = as.numeric(str_sub(step, start = 5, end = -1L))) 
print(dim(dtw_traj_long_rmdup)) # dim = (15440,3)
dtw_traj_sort = dtw_traj_long_rmdup %>% group_by(traj_id) %>% 
  mutate(tot_step = n(), max_step = max(tot_step)) %>% ungroup() %>% 
  select(traj_id, max_step) %>% unique() %>% arrange(desc(max_step)) 
dtw_traj_sort_id = dtw_traj_sort %>% pull(traj_id) 
length(dtw_traj_sort_id) # N = 3272

# COMMAND ----------

# run DTW
install.packages("dtw")
library(dtw)

# COMMAND ----------

# example codes
norm_dist_full = c()
for (i in 1:350) {
  print(i)
  traj_id_ref = dtw_traj_sort_id[i]
  traj_ref = dtw_traj_long_rmdup %>% 
    filter(traj_id == traj_id_ref) %>% pull(ICD_3digit)
  norm_dist_vec = c(rep(NA,i-1))
  for (j in i:length(dtw_traj_sort_id)) {
    traj_i = dtw_traj_sort_id[j]
    traj = dtw_traj_long_rmdup %>% filter(traj_id == traj_i) %>% pull(ICD_3digit)
    step_gap = max(abs(length(traj_ref) - length(traj))-1,1)
    traj_dist_matrix = embed.dist[traj,traj_ref]
    # Perform alignment
    alignment_mvm = dtw(traj_dist_matrix, keep = FALSE, open.end = F, step = asymmetricP05, open.begin = T)
    norm_dist = alignment_mvm$normalizedDistance
    norm_dist_vec = c(norm_dist_vec, norm_dist)
  }
  norm_dist_full = c(norm_dist_full, norm_dist_vec)
  save(norm_dist_full, file = paste0(raw_data_path, "dtw/tmp_results/norm_dist_full_1.rda"))
}

# COMMAND ----------

norm_dist_matrix = matrix(norm_dist_full, nrow = 3272, ncol = 1442)
rownames(norm_dist_matrix) = dtw_traj_sort_id
colnames(norm_dist_matrix) = dtw_traj_sort_id[1831:length(dtw_traj_sort_id)]
save(norm_dist_matrix, file = paste0(raw_data_path, "dtw/tmp_results/norm_dist_matrix_1.rda"))

# COMMAND ----------

# MAGIC %md
# MAGIC ## Combine results together

# COMMAND ----------

for (i in 1:5) {
  load(file = paste0(raw_data_path, "dtw/tmp_results/norm_dist_matrix_", i, ".rda"))
  if (i == 1) {
    norm_dist_matrix_full = norm_dist_matrix
  } else {
    norm_dist_matrix_full = cbind(norm_dist_matrix_full, norm_dist_matrix)
  }
}
print(dim(norm_dist_matrix_full))

# COMMAND ----------

norm_dist_matrix = norm_dist_matrix_full
# print(rownames(norm_dist_matrix) == dtw_traj_sort_id)
# print(colnames(norm_dist_matrix) == dtw_traj_sort_id)
norm_dist_matrix[upper.tri(norm_dist_matrix)] = t(norm_dist_matrix)[upper.tri(norm_dist_matrix)]
print(dim(norm_dist_matrix))
save(norm_dist_matrix, file = paste0(raw_data_path, "dtw/norm_dist_matrix_0313.rda"))

# COMMAND ----------

