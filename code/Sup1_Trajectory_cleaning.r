# Databricks notebook source
# install.packages("igraph")
# library(igraph)
install.packages("R.utils")
library(R.utils)
library(tidyverse)
'%!in%' = function(x,y)!{'%in%'(x,y)}
is_subset = function(v, lst) {
  for (item in lst) {
    if (!identical(v, item) && all(v %in% item)) { return(TRUE) }
  }
  FALSE
}

# COMMAND ----------

raw_data_path = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/mingzhoufu@mednet.ucla.edu/AD_trajectory/data/"
load(file = paste0(raw_data_path, "FG_test/risk_pair_sig_0610.rda"))
print(dim(risk_pair_sig)) # dim = (32670,8)
load(file = paste0(raw_data_path, "patient_data/mod/AD_encounters_filter_0610.rda"))
print(dim(AD_encounters_filter)) # dim = (974718,5)
print(length(unique(AD_encounters_filter$PatientID))) # 24473
head(AD_encounters_filter)

# COMMAND ----------

# MAGIC %md
# MAGIC # 0. Check raw steps

# COMMAND ----------

load(file = paste0(raw_data_path, "FG_test/other_risk/diag_beforeAD_filter_0610.rda"))
print(dim(diag_beforeAD_filter)) # dim = (393601,4)

# COMMAND ----------

# MAGIC %md
# MAGIC ## Move diagnoses at the same time of AD as one step before

# COMMAND ----------

# re-label steps in encounter variables
step_encounters_AD_raw = diag_beforeAD_filter %>% 
  mutate(EncounterDate = as.Date(EncounterDate)) %>% group_by(PatientID) %>% arrange(PatientID, EncounterDate) %>%
  mutate(step = dense_rank(EncounterDate), max_step_raw = max(step)) %>% 
  mutate(mod_step = case_when(
    step == max_step_raw & ICD_3digit != "G30" ~ step - 0.5,
    TRUE ~ step
  )) %>% mutate(step_mod = dense_rank(mod_step), max_step_raw = max(step_mod)) %>% 
  ungroup() %>% unique() %>% select(-c(mod_step, step)) %>% rename("step" = "step_mod")
print(dim(step_encounters_AD_raw)) # dim = (393601,6)

# COMMAND ----------

# check number of steps
step_per_pat_raw = step_encounters_AD_raw %>% select(PatientID, max_step_raw) %>% unique() %>%
  group_by(max_step_raw) %>% summarize(n_patients = n()) %>% arrange(max_step_raw)
print(head(step_per_pat_raw))
print(summary(step_per_pat_raw$max_step_raw))

# COMMAND ----------

# MAGIC %md
# MAGIC Drop N = 1336 people have G30 (AD) diagnosis only.

# COMMAND ----------

# MAGIC %md
# MAGIC ## Preprocess of trajectories

# COMMAND ----------

patient_id_lst_full = step_encounters_AD_raw %>% filter(max_step_raw >= 2) %>% pull(PatientID) %>% unique()
print(length(patient_id_lst_full)) # 23137

# COMMAND ----------

# record patient ids that don't have path to G30 after initial cleaning
exclude_pat_id = c()

# COMMAND ----------

for (i in 1:length(patient_id_lst_full)) {
  if (i %% 1000 == 1) { print(i) }
  patient_id = patient_id_lst_full[i]
  sample_traj = step_encounters_AD_raw %>% filter(PatientID == patient_id)
  # remove ICD codes that are not in the risk pair
  sample_icd = sample_traj$ICD_3digit %>% unique() 
  include_icd = risk_pair_sig %>% filter(exposure_ICD %in% sample_icd & outcome_ICD %in% sample_icd) 
  exclude_icd = setdiff(sample_icd, c(include_icd$exposure_ICD, include_icd$outcome_ICD))
  if ("G30" %in% exclude_icd) {
    exclude_pat_id = c(exclude_pat_id, patient_id)
  } else {
    sample_traj = sample_traj %>% filter(ICD_3digit %!in% exclude_icd) %>% 
      mutate(step = dense_rank(step), max_step_1 = max(step, na.rm = T)) %>% unique()
    if (i == 1) {
      pat_traj_processed = sample_traj
    } else {
      pat_traj_processed = rbind(pat_traj_processed, sample_traj)
    }
  }
}
print(dim(pat_traj_processed)) # dim = (340730,7)
print(length(exclude_pat_id)) # 5585
save(pat_traj_processed, file = paste0(raw_data_path, "patient_data/final/pat_traj_processed_0610.rda"))

# COMMAND ----------

# MAGIC %md
# MAGIC # 1. Deal with people have fewer steps (n = 2)

# COMMAND ----------

load(file = paste0(raw_data_path, "patient_data/final/pat_traj_processed_0610.rda"))
length(unique(pat_traj_processed$PatientID)) # 17552

# COMMAND ----------

few_step_pat_id = pat_traj_processed %>% filter(max_step_1 == 2) %>% pull(PatientID) %>% unique()
print(length(few_step_pat_id)) # 3466
max_step_check = pat_traj_processed %>% filter(max_step_1 == 2) %>% group_by(PatientID) %>%
  summarise(n_diag = n_distinct(ICD_3digit))
table(max_step_check$n_diag)

# COMMAND ----------

# MAGIC %md
# MAGIC ## Short trajectories after expanding - direct impute to multi-steps

# COMMAND ----------

# record patient ids that have too long trajectories (>10) in short-step patients
short_exclude_pat_id = c()
# expand their steps to multiple (based on risk pairs) - strong assumption made (also ran in separate to save time)
for (i in 1:length(few_step_pat_id)) {
  print("===========")
  print(paste0("Patient ", i))
  patient_id = few_step_pat_id[i]
  pat_traj = pat_traj_processed %>% filter(PatientID == few_step_pat_id[i])
  # print(pat_traj)
  node_co_occur = pat_traj %>% arrange(desc(step)) %>% rowid_to_column("id") %>% 
    select(id, ICD_3digit) %>% rename("label" = "ICD_3digit")
  print(paste0("Nodes: ", dim(node_co_occur)[1]))
  if (dim(node_co_occur)[1] <= 10) {
    edges_co_occur = risk_pair_sig %>% 
      filter(exposure_ICD %in% node_co_occur$label & outcome_ICD %in% node_co_occur$label) %>%
      select(exposure_ICD, outcome_ICD) %>% filter(exposure_ICD != "G30") %>%
      left_join(node_co_occur, by = c("outcome_ICD" = "label")) %>% rename(end = outcome_ICD) %>%
      left_join(node_co_occur, by = c("exposure_ICD" = "label")) %>% rename(start = exposure_ICD) %>% select(start, end)
    if (dim(edges_co_occur)[1] > 0) {
      # Create a graph
      g = graph_from_data_frame(edges_co_occur, directed = TRUE)
      # Get all paths from all vertices to vertex F
      allPathsTo_G30 = lapply(1:length(V(g)), function(v) {
        if(V(g)$name[v] != "G30") { # Ensure not to find paths from G30 to G30
          paths = all_simple_paths(g, from = V(g)$name[v], to = "G30")
          return(lapply(paths, function(p) V(g)[p]$name))
        }
      })
      allPathsTo_G30 = Filter(function(x) !is.null(x) && length(x) > 0, allPathsTo_G30) %>% unlist(recursive = FALSE)
      # Filter out elements with length <= 2
      allPathsTo_G30 = Filter(function(x) length(x) >= 2, allPathsTo_G30)
      #===== All possible paths =====#
      full_PathsTo_G30 = allPathsTo_G30
      # Find the maximum length of the paths
      max_length_full = max(sapply(full_PathsTo_G30, length))
      # Pad the shorter paths with NAs to make all paths equal in length
      full_PathsTo_G30 = lapply(full_PathsTo_G30, function(x) {
        length(x) = max_length_full  
        # This pads the vector with NAs if its length is less than max_length
        return(x)
      })
      # Convert the list of paths to a data frame
      paths_pat_full = as.data.frame(do.call(rbind, full_PathsTo_G30), stringsAsFactors = FALSE)
      paths_pat_full$PatientID = patient_id
      paths_pat_full_final = paths_pat_full %>% as.data.frame() %>% unique() 
      # save
      save(paths_pat_full_final, file = paste0(raw_data_path, "trajectory_clean/short_impute/pat_full_new/", i, "_", patient_id, ".rda"))
    } else {
      print("No graph built!")
      short_exclude_pat_id = c(short_exclude_pat_id, patient_id)
    }
  } else {
    print("Too long traj, skip!")
    short_exclude_pat_id = c(short_exclude_pat_id, patient_id)
  }
  save(short_exclude_pat_id, file = paste0(raw_data_path, "trajectory_clean/short_impute/short_exclude_pat_id.rda"))
}

# COMMAND ----------

# combine full results
folder_path = paste0(raw_data_path, "trajectory_clean/short_impute/pat_full_new/") 
# List all .rda files from the folder
files = list.files(path = folder_path, pattern = "\\.rda$")
print(length(files)) # 2790
for (i in 1:length(files)){
  print(i)
  load(paste0(folder_path, files[i]))
  # print(dim(paths_pat_full_final))
  if (i == 1) {
    short_impute = paths_pat_full_final
  } else {
    short_impute = bind_rows(short_impute, paths_pat_full_final) %>% unique()
  }
}
print(dim(short_impute)) # dim = (341868,11)
save(short_impute, file = paste0(raw_data_path, "trajectory_clean/short_impute/short_impute_0610.rda"))

# COMMAND ----------

# MAGIC %md
# MAGIC ## Long-trajectories after expanding - runout time restriction

# COMMAND ----------

load(file = paste0(raw_data_path, "trajectory_clean/short_impute/short_exclude_pat_id.rda"))
print(length(short_exclude_pat_id)) # 676

# COMMAND ----------

# MAGIC %md
# MAGIC Ran this in the background to save time. See details at `~/sep_run/short_long_1`.

# COMMAND ----------

# combine full results
folder_path = paste0(raw_data_path, "trajectory_clean/short_impute/add_long/pat_full/") 
# List all .rda files from the folder
files = list.files(path = folder_path, pattern = "\\.rda$")
print(length(files)) # 360
for (i in 1:length(files)){
  print(i)
  load(paste0(folder_path, files[i]))
  # print(dim(paths_pat_full_final))
  if (i == 1) {
    short_impute_add = paths_pat_full_final
  } else {
    short_impute_add = bind_rows(short_impute_add, paths_pat_full_final) %>% unique()
  }
  tryCatch({
  withTimeout({
    # Place the code that you want to execute within each iteration here
    short_impute_add = short_impute_add %>% filter(!is.na(V3) & is.na(V4)) %>%
      select(PatientID, V1, V2, V3)
    }, timeout = 600)  # Timeout in seconds (600 seconds = 10 minutes)
  }, TimeoutException = function(ex) {
    cat(sprintf("Timeout at iteration %d; moving to next iteration.\n", i))
      
  }, error = function(e) {
    cat("Error in iteration", i, ": ", e$message, "\n") 
  })
  save(short_impute_add, file = paste0(raw_data_path, "trajectory_clean/short_impute/short_impute_add_0610.rda"))
}
print(dim(short_impute_add)) # dim = (2965,4)

# COMMAND ----------

# MAGIC %md
# MAGIC ## Combine results for the short-step people

# COMMAND ----------

load(file = paste0(raw_data_path, "trajectory_clean/short_impute/short_impute_0610.rda"))
print(dim(short_impute)) # dim = (341868,11)
load(file = paste0(raw_data_path, "trajectory_clean/short_impute/short_impute_add_0610.rda"))
print(dim(short_impute_add)) # dim = (2965,4)
short_impute_final = bind_rows(short_impute, short_impute_add) %>% unique()
print(dim(short_impute_final)) # dim = (344833,11)
save(short_impute_final, file = paste0(raw_data_path, "trajectory_clean/short_impute/final/short_impute_final_0610.rda"))

# COMMAND ----------

short_impute_wide_3step = short_impute_final %>% filter(!is.na(V3) & is.na(V4)) %>%
  select(PatientID, V1, V2, V3)
names(short_impute_wide_3step) = c("PatientID", paste0("step", c(1:3)))
print(dim(short_impute_wide_3step)) # dim = (6548,4)
print(length(unique(short_impute_wide_3step$PatientID))) # 1541
save(short_impute_wide_3step, file = paste0(raw_data_path, "trajectory_clean/short_impute/final/short_impute_wide_3step_0610.rda"))
head(short_impute_wide_3step)

# COMMAND ----------

# number of trajectories per person
ttt = short_impute_wide_3step %>% group_by(PatientID) %>% summarise(n = n())
table(ttt$n)

# COMMAND ----------

# number of patients per trajectory
ttt = short_impute_wide_3step %>% group_by(step1, step2, step3) %>% summarise(n_pat = n()) %>% arrange(desc(n_pat))
head(ttt)

# COMMAND ----------

# MAGIC %md
# MAGIC # 2. People with >=3 step trajectories

# COMMAND ----------

full_step_pat_id = pat_traj_processed %>% filter(max_step_1 >= 3) %>% pull(PatientID) %>% unique()
print(length(full_step_pat_id)) # 14086

# COMMAND ----------

# MAGIC %md
# MAGIC ### Identify potential long trajectory patients

# COMMAND ----------

long_traj_pat_id = c()
skip_pat_id = c()

# COMMAND ----------

# Iterate to print potential patients with long trajectories (#step>10)
for (i in 1:length(full_step_pat_id)) {
  if (i %% 500 == 1) { print(i) }
  patient_id = full_step_pat_id[i]
  sample_traj = pat_traj_processed %>% filter(PatientID == patient_id)
  if (dim(sample_traj)[1] > 0) {
    max_step_pat = max(sample_traj$step)
    AD_last_step = sample_traj %>% filter(step == max_step_1 - 1) %>% 
      select(ICD_3digit) %>% rename("last_step" = "ICD_3digit")
    # iterations to clean other steps
    for (s in 2:(max_step_pat-1)) {
      # print(paste0("Step to the last: ", s))
      test_codes_pre = sample_traj %>% filter(step == max_step_1 - s) %>% select(ICD_3digit)
      test_codes = test_codes_pre %>% cross_join(AD_last_step) %>% 
        left_join(risk_pair_sig, by = c("last_step" = "outcome_ICD", "ICD_3digit" = "exposure_ICD")) %>% 
        mutate(step_keep = case_when(
          is.na(coef) ~ "no",
          TRUE ~ "yes"
        )) %>% filter(step_keep == "yes") %>% unique()
      if (nrow(test_codes) > 0) {
        test_codes = test_codes %>% select(ICD_3digit) %>% rename("last_step" = "ICD_3digit") %>% unique()
        AD_last_step = rbind(AD_last_step, test_codes) %>% unique()
      } else {
        sample_traj = sample_traj %>% filter(ICD_3digit %!in% test_codes_pre$ICD_3digit)
      }
    }
    sample_traj_clean = sample_traj %>% group_by(PatientID) %>% 
      mutate(step = dense_rank(step), max_step_1 = max(step)) %>% ungroup() %>% unique()
    # expand to multiple trajectories
    patient_diag_lst = list()
    for (n in 1:max(sample_traj_clean$step)) {
      step_ICD = sample_traj_clean %>% filter(step == max(sample_traj_clean$step)-n+1) %>% pull(ICD_3digit)
      if (n >= 2) {
        unique_elements = unique(unlist(patient_diag_lst))
        for (icd in step_ICD) {
          risk_pair_icd = risk_pair_sig %>% filter(exposure_ICD == icd & outcome_ICD %in% unique_elements) 
          if (nrow(risk_pair_icd) == 0) {
            step_ICD = setdiff(step_ICD, icd)
          }
        }
      }
      patient_diag_lst[[max(sample_traj_clean$step)-n+1]] = step_ICD
    }
    if (length(patient_diag_lst) > 10) {
      long_traj_pat_id = c(long_traj_pat_id, patient_id) %>% unique()
    }
  } else {
    skip_pat_id = c(skip_pat_id, patient_id) %>% unique()
  }
}

# COMMAND ----------

print(length(long_traj_pat_id)) # 5278
print(length(skip_pat_id)) # 0
save(long_traj_pat_id, file = paste0(raw_data_path, "trajectory_clean/no_impute/long_traj_pat_id.rda"))
save(skip_pat_id, file = paste0(raw_data_path, "trajectory_clean/no_impute/skip_pat_id.rda"))

# COMMAND ----------

load(file = paste0(raw_data_path, "trajectory_clean/no_impute/long_traj_pat_id.rda"))
load(file = paste0(raw_data_path, "trajectory_clean/no_impute/skip_pat_id.rda"))
print(length(long_traj_pat_id)) # 5278
print(length(skip_pat_id)) # 0
load(file = paste0(raw_data_path, "patient_data/final/pat_traj_processed_0610.rda"))

# COMMAND ----------

# MAGIC %md
# MAGIC ### Long trajectories: find potential short add back

# COMMAND ----------

short_add_pat_id = c()
sep_run_pat_id = c()
long_remove_pat_id = c()

# COMMAND ----------

load(file = paste0(raw_data_path, "trajectory_clean/no_impute/long_traj/sep_run_pat_id.rda"))
load(file = paste0(raw_data_path, "trajectory_clean/no_impute/long_traj/short_add_pat_id.rda"))
load(file = paste0(raw_data_path, "trajectory_clean/no_impute/long_traj/long_remove_pat_id.rda"))
index_sep = which(long_traj_pat_id %in% sep_run_pat_id) %>% max()
index_short = which(long_traj_pat_id %in% short_add_pat_id) %>% max()
index_long = which(long_traj_pat_id %in% long_remove_pat_id) %>% max()
max_index = max(index_sep, index_short, index_long)
error_index = max_index + 1
print(error_index)

# COMMAND ----------

# Iterate to print potential patients with long trajectories (#step>10)
for (i in (error_index+1):length(long_traj_pat_id)) {
  print(i)
  patient_id = long_traj_pat_id[i]
  sample_traj = pat_traj_processed %>% filter(PatientID == patient_id)
  if (dim(sample_traj)[1] > 0) {
    max_step_pat = max(sample_traj$step)
    AD_last_step = sample_traj %>% filter(step == max_step_1 - 1) %>% 
      select(ICD_3digit) %>% rename("last_step" = "ICD_3digit")
    # iterations to clean other steps
    for (s in 2:(max_step_pat-1)) {
      # print(paste0("Step to the last: ", s))
      test_codes_pre = sample_traj %>% filter(step == max_step_1 - s) %>% select(ICD_3digit)
      test_codes = test_codes_pre %>% cross_join(AD_last_step) %>% 
        left_join(risk_pair_sig, by = c("last_step" = "outcome_ICD", "ICD_3digit" = "exposure_ICD")) %>% 
        mutate(step_keep = case_when(
          is.na(coef) ~ "no",
          TRUE ~ "yes"
        )) %>% filter(step_keep == "yes") %>% unique()
      if (nrow(test_codes) > 0) {
        test_codes = test_codes %>% select(ICD_3digit) %>% rename("last_step" = "ICD_3digit") %>% unique()
        AD_last_step = rbind(AD_last_step, test_codes) %>% unique()
      } else {
        sample_traj = sample_traj %>% filter(ICD_3digit %!in% test_codes_pre$ICD_3digit)
      }
    }
    sample_traj_clean = sample_traj %>% group_by(PatientID) %>% 
      mutate(step = dense_rank(step), max_step_1 = max(step)) %>% ungroup() %>% unique()
    # expand to multiple trajectories
    patient_diag_lst = list()
    for (n in 1:max(sample_traj_clean$step)) {
      step_ICD = sample_traj_clean %>% filter(step == max(sample_traj_clean$step)-n+1) %>% pull(ICD_3digit)
      if (n >= 2) {
        unique_elements = unique(unlist(patient_diag_lst))
        for (icd in step_ICD) {
          risk_pair_icd = risk_pair_sig %>% filter(exposure_ICD == icd & outcome_ICD %in% unique_elements) 
          if (nrow(risk_pair_icd) == 0) {
            step_ICD = setdiff(step_ICD, icd)
          }
        }
      }
      patient_diag_lst[[max(sample_traj_clean$step)-n+1]] = step_ICD
    }
    patient_diag_lst = patient_diag_lst[sapply(patient_diag_lst, length) > 0]
    tryCatch({
      withTimeout({
        expanded_list = expand.grid(patient_diag_lst) %>% t() %>% as.data.frame() %>%
          mutate(step = row_number(), PatientID = patient_id)
        # If the length >10,000 and <50,000, we save the patients' matrix and run them separately, recorded their ids as sep_run
        if (length(expanded_list) > 10000 & length(expanded_list) < 50000) {
          sep_run_pat_id = c(sep_run_pat_id, patient_id)
          save(sep_run_pat_id, file = paste0(raw_data_path, "trajectory_clean/no_impute/long_traj/sep_run_pat_id.rda"))
          save(expanded_list, file = paste0(raw_data_path, "trajectory_clean/no_impute/long_traj/tmp_expand/", patient_id, ".rda"))
        } else if (length(expanded_list) <= 10000) {
          short_add_pat_id = c(short_add_pat_id, patient_id)
          save(short_add_pat_id, file = paste0(raw_data_path, "trajectory_clean/no_impute/long_traj/short_add_pat_id.rda"))
        } else {
          # If the length >100,000 or it took a long time to expand, we record their ids as long_traj (those will be skipped for cleaning)
          long_remove_pat_id = c(long_remove_pat_id, patient_id)
          save(long_remove_pat_id, file = paste0(raw_data_path, "trajectory_clean/no_impute/long_traj/long_remove_pat_id.rda"))
        }
      }, timeout = 600)
    }, TimeoutException = function(ex) {
      cat(sprintf("Timeout at iteration %d; moving to next iteration.\n", i))
      long_remove_pat_id = c(long_remove_pat_id, patient_id)
      save(long_remove_pat_id, file = paste0(raw_data_path, "trajectory_clean/no_impute/long_traj/long_remove_pat_id.rda"))
    }, error = function(e) {
      cat("Error in iteration", i, ": ", e$message, "\n")
      long_remove_pat_id = c(long_remove_pat_id, patient_id)
      save(long_remove_pat_id, file = paste0(raw_data_path, "trajectory_clean/no_impute/long_traj/long_remove_pat_id.rda"))
    })
  } 
}

# COMMAND ----------

load(file = paste0(raw_data_path, "trajectory_clean/no_impute/long_traj/sep_run_pat_id.rda"))
print(length(sep_run_pat_id))
load(file = paste0(raw_data_path, "trajectory_clean/no_impute/long_traj/long_remove_pat_id.rda"))
print(length(long_remove_pat_id))
load(file = paste0(raw_data_path, "trajectory_clean/no_impute/long_traj/short_add_pat_id.rda"))
print(length(short_add_pat_id))

# COMMAND ----------

# MAGIC %md
# MAGIC ### Short trajectories
# MAGIC We just ran code directly

# COMMAND ----------

load(file = paste0(raw_data_path, "patient_data/final/pat_traj_processed_0610.rda"))
full_step_pat_id = pat_traj_processed %>% filter(max_step_1 >= 3) %>% pull(PatientID) %>% unique()
print(length(full_step_pat_id)) # 14086
load(file = paste0(raw_data_path, "trajectory_clean/no_impute/long_traj/short_add_pat_id.rda"))
print(length(short_add_pat_id)) # 4776
load(file = paste0(raw_data_path, "trajectory_clean/no_impute/long_traj_pat_id.rda"))
print(length(long_traj_pat_id)) # 5278
risk_AD_sig = risk_pair_sig %>% filter(outcome_ICD == "G30") 
AD_risk_ICD = risk_AD_sig %>% pull(exposure_ICD) %>% unique() # N = 16

# COMMAND ----------

short_traj_id = c(setdiff(full_step_pat_id, long_traj_pat_id), short_add_pat_id)
print(length(short_traj_id)) # 13584

# COMMAND ----------

# MAGIC %md
# MAGIC We ran this sepearately to save time, see details at `~/short_run_1k`

# COMMAND ----------

# combine full results
folder_path = paste0(raw_data_path, "trajectory_clean/no_impute/short_traj/pat_filter/") 
# List all .rda files from the folder
files = list.files(path = folder_path, pattern = "\\.rda$")
print(length(files)) # 4191
num_files = sapply(strsplit(files, "_"), `[`, 1) %>% as.numeric() %>% sort()

# COMMAND ----------

for (i in 1:length(files)){
  if (i %% 100 == 1) { print(i) }
  load(paste0(folder_path, files[i]))
  # print(dim(paths_pat_filter_final))
  if (i == 1) {
    short_full = paths_pat_filter_final
  } else {
    short_full = bind_rows(short_full, paths_pat_filter_final) %>% unique()
  }
  tryCatch({
  withTimeout({
    # Place the code that you want to execute within each iteration here
    short_full = short_full %>% filter(!is.na(V3)) 
    }, timeout = 600)  # Timeout in seconds (600 seconds = 10 minutes)
  }, TimeoutException = function(ex) {
    cat(sprintf("Timeout at iteration %d; moving to next iteration.\n", i))
      
  }, error = function(e) {
    cat("Error in iteration", i, ": ", e$message, "\n") 
  })
  save(short_full, file = paste0(raw_data_path, "trajectory_clean/no_impute/short_traj/short_full_0622.rda"))
}
print(dim(short_full)) # dim = (79433,18)

# COMMAND ----------

load(file = paste0(raw_data_path, "trajectory_clean/no_impute/short_traj/short_full_0622.rda"))
tail(short_full)

# COMMAND ----------

short_full %>% filter(!is.na(V3)) %>% pull(PatientID) %>% unique() %>% length()

# COMMAND ----------

# MAGIC %md
# MAGIC #### Expand 2-step-after-cleaning patients

# COMMAND ----------

# MAGIC %md
# MAGIC We ran this sepearately to save time, see details at `~/short_to_impute`

# COMMAND ----------

load(file = paste0(raw_data_path, "trajectory_clean/short_impute/final/short_impute_wide_3step_0610.rda"))
print(dim(short_impute_wide_3step)) # dim = (6548,4)
head(short_impute_wide_3step)

# COMMAND ----------

load(file = paste0(raw_data_path, "trajectory_clean/no_impute/short_impute/long_impute_wide_3step_0622.rda"))
impute_wide_3step = rbind(short_impute_wide_3step, long_impute_wide_3step)
print(dim(impute_wide_3step))# dim = (6591,4)
print(length(unique(impute_wide_3step$PatientID))) # 1571
save(impute_wide_3step, file = paste0(raw_data_path, "trajectory_clean/final/impute_wide_3step_0622.rda"))

# COMMAND ----------

load(file = paste0(raw_data_path, "trajectory_clean/no_impute/short_traj/short_full_0622.rda"))
print(dim(short_full)) # dim = (79433,18)
print(length(unique(short_full$PatientID))) # 4191
head(short_full)

# COMMAND ----------

paths_pat_combine = short_full %>% filter(!is.na(V3)) %>% unique() %>% 
  select(PatientID, V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17)
names(paths_pat_combine) = c("PatientID", paste0("step", c(1:17)))
patient_traj_wide = bind_rows(paths_pat_combine, impute_wide_3step)
print(dim(patient_traj_wide)) # dim = (86024,18)
print(length(unique(patient_traj_wide$PatientID)))# 5762
save(patient_traj_wide, file = paste0(raw_data_path, "trajectory_clean/final/patient_traj_wide_0622.rda"))

# COMMAND ----------

head(patient_traj_wide)

# COMMAND ----------

