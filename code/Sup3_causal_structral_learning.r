# Databricks notebook source
install.packages("cpp11")
install.packages("/Workspace/Users/mingzhoufu@mednet.ucla.edu/R-lib/igraph_1.6.0.tar.gz", repos = NULL, type = "source")
install.packages("/Workspace/Users/mingzhoufu@mednet.ucla.edu/R-lib/BiocGenerics_0.48.1.tar.gz", repos = NULL, type = "source")
install.packages("/Workspace/Users/mingzhoufu@mednet.ucla.edu/R-lib/graph_1.80.0.tar.gz", repos = NULL, type = "source")
install.packages("/Workspace/Users/mingzhoufu@mednet.ucla.edu/R-lib/BH_1.84.0-0.tar.gz", repos = NULL, type = "source")
install.packages("/Workspace/Users/mingzhoufu@mednet.ucla.edu/R-lib/RBGL_1.78.0.tar.gz", repos = NULL, type = "source")
install.packages("/Workspace/Users/mingzhoufu@mednet.ucla.edu/R-lib/BiocManager_1.30.22.tar.gz", repos = NULL, type = "source")
install.packages("/Workspace/Users/mingzhoufu@mednet.ucla.edu/R-lib/ggm_2.5.1.tar.gz", repos = NULL, type = "source")
install.packages("pcalg")
library(pcalg)
library(tidyverse)
raw_data_path = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/AD_trajectory/data/"

# COMMAND ----------

# MAGIC %md
# MAGIC # GES causal structural learning run separately

# COMMAND ----------

# MAGIC %md
# MAGIC ## AD causal links

# COMMAND ----------

load(file = paste0(raw_data_path, "causal/AD_links/combined_datasets_case_control.rda"))

# COMMAND ----------

i = 1 # 1-10
test_df = combined_datasets[[i]]
# learn causal structure with GES
score = new("GaussL0penObsScore", data = as.matrix(test_df))
ges.fit = ges(score)
save(ges.fit, file = paste0(raw_data_path, "causal/AD_links/ges.fit_", i, ".rda"))

# COMMAND ----------

# MAGIC %md
# MAGIC ## Causal links other than AD

# COMMAND ----------

load(file = paste0(raw_data_path, "causal/other_links/combined_datasets_AD_only.rda"))

# COMMAND ----------

i = 1 # 1-10
test_df = combined_datasets[[i]]
# learn causal structure with GES
score = new("GaussL0penObsScore", data = as.matrix(test_df))
ges.fit = ges(score)
save(ges.fit, file = paste0(raw_data_path, "causal/other_links/ges.fit_", i, ".rda"))

# COMMAND ----------

# MAGIC %md
# MAGIC # Clean GES results

# COMMAND ----------

# MAGIC %md
# MAGIC ## AD causal links

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
print(dim(causal_full)) # dim = (111290,3)
save(causal_full, file = paste0(raw_data_path, "causal/AD_links/causal_full_AD.rda"))

# COMMAND ----------

# MAGIC %md
# MAGIC ## Causal links other than AD

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
print(dim(causal_full)) # dim = (106847,3)
save(causal_full, file = paste0(raw_data_path, "causal/other_links/causal_full_other.rda"))

# COMMAND ----------

