# Databricks notebook source
# MAGIC %md
# MAGIC #### Last updated: 06/22/2024

# COMMAND ----------

install.packages("R.utils")
library(R.utils)
# install.packages("dtw")
# library(dtw)
library(tidyverse)
'%!in%' = function(x,y)!{'%in%'(x,y)}
is_subset = function(v, lst) {
  for (item in lst) {
    if (!identical(v, item) && all(v %in% item)) { return(TRUE) }
  }
  FALSE
}
raw_data_path = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/mingzhoufu@mednet.ucla.edu/AD_trajectory/data/"
load(file = paste0(raw_data_path, "FG_test/risk_pair_sig_0610.rda"))
print(dim(risk_pair_sig)) # dim = (32670,8)

# COMMAND ----------

# MAGIC %md
# MAGIC # 1. AD patients' trajectory cleaning

# COMMAND ----------

# MAGIC %md
# MAGIC See details at `~/Sup1_traj_cleaning`

# COMMAND ----------

load(file = paste0(raw_data_path, "trajectory_clean/final/patient_traj_wide_0622.rda"))
print(dim(patient_traj_wide)) # dim = (86024,18)
print(length(unique(patient_traj_wide$PatientID))) # 5762

# COMMAND ----------

# MAGIC %md
# MAGIC ## Restrict max step of trajectory -- 9 step

# COMMAND ----------

patient_traj_unique = patient_traj_wide %>% select(-c(PatientID)) %>% unique() %>% mutate(TrajID = row_number()) 
print(dim(patient_traj_unique)) # dim = (73948,18)
patient_traj_wide_add = patient_traj_wide %>% left_join(patient_traj_unique)
print(dim(patient_traj_wide_add)) # dim = (86024,19)
head(patient_traj_wide_add)

# COMMAND ----------

# Counting the number of missing values in each column
missing_values = colSums(is.na(patient_traj_unique))
# Printing the number of missing values for each column
print("Number of missing values in each column:")
print(missing_values)

# COMMAND ----------

dtw_traj_long_rmdup = patient_traj_unique %>% 
  pivot_longer(cols = setdiff(names(patient_traj_unique), "TrajID"), names_to = "step", values_to = "ICD_3digit") %>% drop_na() %>% 
  mutate(step = as.numeric(str_sub(step, start = 5, end = -1L))) 
print(dim(dtw_traj_long_rmdup)) # dim = (630577,3)
head(dtw_traj_long_rmdup)

# COMMAND ----------

max_step_traj = dtw_traj_long_rmdup %>% group_by(TrajID) %>% summarise(traj_step = max(step)) 
pat_step_summary = patient_traj_wide_add %>% select(PatientID, TrajID) %>% left_join(max_step_traj) %>% group_by(PatientID) %>%
  summarise(max_step_pat = max(traj_step), min_step_pat = min(traj_step))
  print(summary(pat_step_summary$min_step_pat))
print(table(pat_step_summary$min_step_pat))
dim(pat_step_summary)[1]*0.01 # 99% of all patients have a minimum step no more than 9-step

# COMMAND ----------

traj_wide_step_filter = patient_traj_wide %>% filter(is.na(step10)) %>% select(-c(paste0("step", c(10:17))))
print(dim(traj_wide_step_filter)) # dim = (59801,10)
print(length(unique(traj_wide_step_filter$PatientID)))
head(traj_wide_step_filter) # 5726

# COMMAND ----------

# MAGIC %md
# MAGIC ## Add back in excluded patients (with >9 step only)

# COMMAND ----------

exclude_pat_id = setdiff(unique(patient_traj_wide$PatientID), unique(traj_wide_step_filter$PatientID)) 
print(length(exclude_pat_id)) # N = 36

# COMMAND ----------

# combine full results
folder_path = paste0(raw_data_path, "trajectory_clean/no_impute/short_traj/pat_full/") 
# List all .rda files from the folder
files = list.files(path = folder_path, pattern = "\\.rda$")
print(length(files)) # 5246
pat_id_files_pre = sapply(strsplit(files, "_"), `[`, 2)
# Function to remove the last 4 characters from each element
remove_last_4 = function(x) {
  substr(x, 1, nchar(x) - 4)
}
# Apply the function to the vector
pat_id_files = sapply(pat_id_files_pre, remove_last_4)
select_files_index = which(pat_id_files %in% exclude_pat_id)
print(length(select_files_index)) # 36

# COMMAND ----------

for (i in select_files_index){
  load(paste0(folder_path, files[i]))
  # print(dim(paths_pat_filter_final))
  if (i == select_files_index[1]) {
    short_add = paths_pat_full_final
  } else {
    short_add = bind_rows(short_add, paths_pat_full_final) %>% unique()
  }
  tryCatch({
  withTimeout({
    # Place the code that you want to execute within each iteration here
    short_add = short_add %>% filter(!is.na(V3)) 
    }, timeout = 600)  # Timeout in seconds (600 seconds = 10 minutes)
  }, TimeoutException = function(ex) {
    cat(sprintf("Timeout at iteration %d; moving to next iteration.\n", i))
      
  }, error = function(e) {
    cat("Error in iteration", i, ": ", e$message, "\n") 
  })
}
print(dim(short_add)) # dim = (214135,17)

# COMMAND ----------

short_add_filter = short_add %>% filter(!is.na(V9) & is.na(V10)) %>% select(PatientID, c(paste0("V", c(1:9)))) %>% unique()
names(short_add_filter) = c("PatientID", c(paste0("step", c(1:9))))
print(dim(short_add_filter)) # dim = (45229,10)
head(short_add_filter)

# COMMAND ----------

traj_wide_step_filter_add = rbind(traj_wide_step_filter, short_add_filter) %>% unique()
print(dim(traj_wide_step_filter_add)) # dim = (105030,10)
length(unique(traj_wide_step_filter_add$PatientID)) # N = 5762

# COMMAND ----------

save(traj_wide_step_filter_add, file = paste0(raw_data_path, "trajectory_clean/final/traj_wide_step_filter_add.rda"))

# COMMAND ----------

# MAGIC %md
# MAGIC ## Restrict #traj per patient -- up to 3 (median)

# COMMAND ----------

load(file = paste0(raw_data_path, "trajectory_clean/final/traj_wide_step_filter_add.rda"))
print(dim(traj_wide_step_filter_add)) # dim = (105030,10)
length(unique(traj_wide_step_filter_add$PatientID)) # N = 5762

# COMMAND ----------

traj_unique_raw = traj_wide_step_filter_add %>% select(-PatientID) %>% unique()
print(dim(traj_unique_raw)) # dim = (92954,9)

# COMMAND ----------

traj_per_pat = traj_wide_step_filter_add %>% group_by(PatientID) %>% summarise(n_traj = n())
summary(traj_per_pat$n_traj)

# COMMAND ----------

sum(traj_per_pat$n_traj > 3)/dim(traj_per_pat)[1] # ~50% of patients need to do the sampling of their trajectories

# COMMAND ----------

# MAGIC %md
# MAGIC ## Sampling from the trajectories 

# COMMAND ----------

sample_pat_id = traj_per_pat %>% filter(n_traj > 3) %>% pull(PatientID) %>% unique()
print(length(sample_pat_id)) # 2846
save(sample_pat_id, file = paste0(raw_data_path, "trajectory_clean/final/sample_pat_id.rda"))

# COMMAND ----------

# MAGIC %md
# MAGIC For those patients with extreme large number of trajectories, to get 3 representative trajectories.
# MAGIC See details at `~/traj_sample_1`

# COMMAND ----------

load(file = paste0(raw_data_path, "trajectory_clean/final/sampled_pat_traj_final.rda"))
print(dim(sampled_pat_traj_final))
head(sampled_pat_traj_final)

# COMMAND ----------

# MAGIC %md
# MAGIC ## Finalize patient trajectories

# COMMAND ----------

pat_traj_full = traj_wide_step_filter_add %>% filter(PatientID %!in% sample_pat_id)
pat_traj_final = rbind(pat_traj_full, sampled_pat_traj_final) %>% as.data.frame() %>% unique()
print(dim(pat_traj_final)) # dim = (13409, 10)
save(pat_traj_final, file = paste0(raw_data_path, "trajectory_clean/final/pat_traj_final_0622.rda"))

# COMMAND ----------

traj_per_pat = pat_traj_final %>% group_by(PatientID) %>% summarise(n_traj = n())
print(dim(traj_per_pat)) # dim = (5762, 2)
summary(traj_per_pat$n_traj)

# COMMAND ----------

patient_traj_unique = pat_traj_final %>% select(-PatientID) %>% unique() %>% mutate(traj_id = row_number()) 
print(dim(patient_traj_unique)) # dim = (6794,10)
dtw_traj_long_rmdup = patient_traj_unique %>% 
  pivot_longer(cols = setdiff(names(patient_traj_unique), "traj_id"), names_to = "step", values_to = "ICD_3digit") %>% drop_na() %>% 
  mutate(step = as.numeric(str_sub(step, start = 5, end = -1L))) 
print(dim(dtw_traj_long_rmdup)) # dim = (40301,3)

# COMMAND ----------

step_per_traj = dtw_traj_long_rmdup %>% group_by(traj_id) %>%
  summarise(max_step = max(step))
summary(step_per_traj$max_step)

# COMMAND ----------

# MAGIC %md
# MAGIC # 2. DTW alignment

# COMMAND ----------

# MAGIC %md
# MAGIC ## Data preparation

# COMMAND ----------

patient_traj_unique = pat_traj_final %>% select(-PatientID) %>% unique() %>% mutate(traj_id = row_number()) 
print(dim(patient_traj_unique)) # dim = (6794,10)

# COMMAND ----------

dtw_traj_long_rmdup = patient_traj_unique %>% 
  pivot_longer(cols = setdiff(names(patient_traj_unique), "traj_id"), names_to = "step", values_to = "ICD_3digit") %>% drop_na() %>% 
  mutate(step = as.numeric(str_sub(step, start = 5, end = -1L))) 
print(dim(dtw_traj_long_rmdup)) # dim = (40301,3)
dtw_traj_sort = dtw_traj_long_rmdup %>% group_by(traj_id) %>% 
  mutate(tot_step = n(), max_step = max(tot_step)) %>% ungroup() %>% 
  select(traj_id, max_step) %>% unique() %>% arrange(desc(max_step)) 
dtw_traj_sort_id = dtw_traj_sort %>% pull(traj_id) 
length(dtw_traj_sort_id) # N = 6794

# COMMAND ----------

dtw_traj_sort = dtw_traj_long_rmdup %>% group_by(traj_id) %>% 
  mutate(tot_step = n(), max_step = max(tot_step)) %>% ungroup() %>% 
  select(traj_id, max_step) %>% unique() %>% arrange(desc(max_step)) 
dtw_traj_sort_id = dtw_traj_sort %>% pull(traj_id) 
length(dtw_traj_sort_id) # N = 6794

# COMMAND ----------

save(dtw_traj_sort_id, file = paste0(raw_data_path, "trajectory_clean/final/dtw_traj_sort_id_0622.rda"))
save(dtw_traj_long_rmdup, file = paste0(raw_data_path, "trajectory_clean/final/dtw_traj_long_rmdup_0622.rda"))

# COMMAND ----------

# MAGIC %md
# MAGIC ## DTW distance matrix
# MAGIC We also ran this separately, see details at `~/DTW_bg_1`.

# COMMAND ----------

load(file = paste0(raw_data_path, "dtw/norm_dist_matrix_0622.rda"))
print(dim(norm_dist_matrix)) # dim = (6794,6794)

# COMMAND ----------

# MAGIC %md
# MAGIC # 3. Clustering

# COMMAND ----------

# MAGIC %md
# MAGIC ## Determine optimal number of clusters

# COMMAND ----------

# MAGIC %md
# MAGIC We also ran this separately, see details at `~/cluster_hc`.

# COMMAND ----------

install.packages("clusterSim")
library(clusterSim)
install.packages("clusterCrit")
library(clusterCrit)
install.packages("cluster")
library(cluster)
install.packages("mclust")
library(mclust)
library(tidyverse)

# COMMAND ----------

# load in distance matrix
raw_data_path = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/mingzhoufu@mednet.ucla.edu/AD_trajectory/data/"
load(file = paste0(raw_data_path, "dtw/norm_dist_matrix_0622.rda"))
dim(norm_dist_matrix) # dim = (6794,6794)

# COMMAND ----------

load(file = paste0(raw_data_path, "clusters/tmp/hc_result.rda"))
load(file = paste0(raw_data_path, "clusters/tmp/pam_result.rda"))
load(file = paste0(raw_data_path, "clusters/tmp/knn_result.rda"))

# COMMAND ----------

hc_result

# COMMAND ----------

pam_result

# COMMAND ----------

knn_result

# COMMAND ----------

# MAGIC %md
# MAGIC According to results from all three methods， we decide to set N_cluster = 4.

# COMMAND ----------

# MAGIC %md
# MAGIC ## Combine clustering results

# COMMAND ----------

# Hierarchical (N = 4 clusters)
final_sample_hc = hclust(as.dist(norm_dist_matrix), method = "ward.D2")
cluster_hc = cutree(final_sample_hc, k = 4)
# PAM (N = 4 clusters)
final_sample_pam = pam(as.dist(norm_dist_matrix), k = 4)
# K-NN (N = 4 clusters)
final_sample_knn = kmeans(norm_dist_matrix, 4)

# COMMAND ----------

# Compare clusters
cluster_comparison = data.frame(
  traj_id = colnames(norm_dist_matrix),
  Kmeans = final_sample_knn$cluster,
  Hierarchical = cluster_hc,
  PAM = final_sample_pam$clustering
)
cluster_comparison$traj_id = as.numeric(cluster_comparison$traj_id)
print(dim(cluster_comparison)) # dim = (7839,4)
head(cluster_comparison)

# COMMAND ----------

# Create confusion matrix for each pair of clustering results
confusion_matrix_kmeans_hc = table(cluster_comparison$Kmeans, cluster_comparison$Hierarchical)
confusion_matrix_kmeans_pam = table(cluster_comparison$Kmeans, cluster_comparison$PAM)
confusion_matrix_hc_pam = table(cluster_comparison$Hierarchical, cluster_comparison$PAM)
# Print confusion matrices
print("Confusion Matrix between K-means and Hierarchical Clustering:")
print(confusion_matrix_kmeans_hc)
print("Confusion Matrix between K-means and PAM Clustering:")
print(confusion_matrix_kmeans_pam)
print("Confusion Matrix between Hierarchical and PAM Clustering:")
print(confusion_matrix_hc_pam)

# COMMAND ----------

# Calculate adjusted Rand index for each pair of clustering results
ari_kmeans_hc = adjustedRandIndex(cluster_comparison$Kmeans, cluster_comparison$Hierarchical)
ari_kmeans_pam = adjustedRandIndex(cluster_comparison$Kmeans, cluster_comparison$PAM)
ari_hc_pam = adjustedRandIndex(cluster_comparison$Hierarchical, cluster_comparison$PAM)

# Print adjusted Rand indices
print(paste("Adjusted Rand Index between K-means and Hierarchical Clustering:", round(ari_kmeans_hc, 3)))
print(paste("Adjusted Rand Index between K-means and PAM Clustering:", round(ari_kmeans_pam, 3)))
print(paste("Adjusted Rand Index between Hierarchical and PAM Clustering:", round(ari_hc_pam, 3)))

# COMMAND ----------

table(cluster_comparison$Kmeans)

# COMMAND ----------

load(file = paste0(raw_data_path, "trajectory_clean/final/pat_traj_final_0622.rda"))
patient_traj_unique = pat_traj_final %>% dplyr::select(-PatientID) %>% unique() %>% mutate(traj_id = row_number()) 
print(dim(patient_traj_unique)) # dim = (6794,10)
# add cluster info
patient_traj_cluster = patient_traj_unique %>% left_join(cluster_comparison)
print(dim(patient_traj_cluster)) # dim = (6794,13)
pat_traj_final_cluster = pat_traj_final %>% left_join(patient_traj_cluster)
print(dim(pat_traj_final)) # dim = (13409,10)
print(dim(pat_traj_final_cluster)) # dim = (13409,14)
save(pat_traj_final_cluster, file = paste0(raw_data_path, "patient_data/final/pat_traj_final_cluster_0622.rda"))
head(pat_traj_final_cluster)

# COMMAND ----------

# cluster summary
cluster_summary = pat_traj_final_cluster %>% dplyr::select(PatientID, traj_id, Kmeans) %>% group_by(Kmeans) %>% 
  summarise(n_pat = n_distinct(PatientID), n_traj = n_distinct(traj_id)) %>% arrange(desc(n_pat))
cluster_summary

# COMMAND ----------

pat_cluster_summary = pat_traj_final_cluster %>% dplyr::select(PatientID, Kmeans) %>% 
  unique() %>% group_by(PatientID) %>% summarise(n_cluster = n()) %>% 
  arrange(desc(n_cluster))
table(pat_cluster_summary$n_cluster)

# COMMAND ----------

