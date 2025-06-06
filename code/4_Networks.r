# Databricks notebook source
# MAGIC %md
# MAGIC #### Last updated: 06/24/2024

# COMMAND ----------

install.packages("igraph")
library(igraph)
install.packages("tidygraph")
library(tidygraph)
install.packages("ggraph")
library(ggraph)
install.packages("backbone")
library(backbone)
install.packages("icd.data")
library(icd.data)
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
print(dim(risk_pair_sig)) # dim = (31703,8)
icd10cm2016_short = icd10cm2016 %>% select(code, long_desc) %>% mutate(code = as.character(code))

# COMMAND ----------

load(file = paste0(raw_data_path, "patient_data/final/pat_traj_final_cluster_0622.rda"))

# COMMAND ----------

# MAGIC %md
# MAGIC # 1. Full network

# COMMAND ----------

# MAGIC %md
# MAGIC ## Raw network
# MAGIC Build raw network of all AD patients (before clustering)

# COMMAND ----------

pat_cluster_i_long = pat_traj_final_cluster %>% 
  select(PatientID, paste0("step", c(1:9))) %>% 
  pivot_longer(cols = paste0("step", c(1:9)), names_to = "step", values_to = "ICD_3digit") %>% drop_na() %>% 
  mutate(step = as.numeric(str_sub(step, start = 5, end = -1L))) 
id_lst = unique(pat_cluster_i_long$PatientID)
print(length(id_lst)) # 5762
for (n in 1:length(id_lst)) {
  patient_id = id_lst[n]
  # print(paste0("Patinet ", n, ": ", patient_id))
  pat_rec = pat_cluster_i_long %>% filter(PatientID == patient_id) %>% 
    group_by(step) %>% mutate(traj = row_number()) %>% ungroup()
  n_traj = max(pat_rec$traj)
  for (i in 1:n_traj) {
    icd_vec = pat_rec %>% filter(traj == i) %>% pull(ICD_3digit)
    start_end_sub = data.frame(start = icd_vec[-length(icd_vec)], end = icd_vec[-1]) %>% 
      mutate(PatientID = patient_id, traj = i, n = 1)
    if (i == 1) {patient_start_end = start_end_sub}
    else {patient_start_end = rbind(patient_start_end, start_end_sub)}
  }
  if (n == 1) {start_end_full = patient_start_end}
  else {start_end_full = rbind(start_end_full, patient_start_end)}
}
pat_start_end_summary = start_end_full %>% select(-traj) %>% unique() %>% 
  group_by(start, end) %>% summarise(n_patient = n()) %>% filter(start != "G30") %>% ungroup()
print(dim(pat_start_end_summary)) # dim = (7930,3)
save(pat_start_end_summary, file = paste0(raw_data_path, "patient_data/final/pat_start_end_summary_full_0624.rda"))

# COMMAND ----------

load(file = paste0(raw_data_path, "patient_data/final/pat_start_end_summary_full_0624.rda"))
pat_start_end_summary %>% filter(n_patient >= 500) %>% arrange(desc(n_patient))

# COMMAND ----------

# filter based on significant pairs and #pat followed (>=0.5%)
pat_start_end_summary_sig = pat_start_end_summary %>% 
  left_join(risk_pair_sig %>% select(exposure_ICD, outcome_ICD, coef) %>% 
  rename(start = exposure_ICD, end = outcome_ICD), by = c("start", "end")) %>% 
  filter(!is.na(coef)) %>% select(-coef) %>% filter(n_patient >= 6984*0.005)
print(dim(pat_start_end_summary_sig)) # dim = (108,3)
# Creating a weighted graph
g = graph_from_data_frame(data.frame(
  from = pat_start_end_summary_sig$start,
  to = pat_start_end_summary_sig$end,
  weight = pat_start_end_summary_sig$n_patient), directed = TRUE)
print(length(V(g)$name)) # 45
# Get all paths from all vertices to vertex G30 (remove edges that cannot lead to G30 b/c we only care about paths to G30)
allPathsTo_G30 = lapply(1:length(V(g)), function(v) {
  if(V(g)$name[v] != "G30") { # Ensure not to find paths from G30 to G30
    paths = all_simple_paths(g, from = V(g)$name[v], to = "G30")
    return(lapply(paths, function(p) V(g)[p]$name))
  }
})
allPathsTo_G30 = Filter(function(x) !is.null(x) && length(x) > 0, allPathsTo_G30)
allPathsTo_G30 = unlist(allPathsTo_G30, recursive = FALSE)
# Find the maximum length of the paths
max_length = max(sapply(allPathsTo_G30, length))
# Pad the shorter paths with NAs to make all paths equal in length
allPathsTo_G30 = lapply(allPathsTo_G30, function(x) {
  length(x) = max_length  
  return(x)
})
# Convert the list of paths to a data frame
paths_df = as.data.frame(do.call(rbind, allPathsTo_G30), stringsAsFactors = FALSE) %>% filter(!is.na(V3))
print(dim(paths_df)) # dim = (13334,13)

# COMMAND ----------

# restrict to nodes to G30
node_to_G30 = unique(c(paths_df$V1, paths_df$V2, paths_df$V3, paths_df$V4, paths_df$V5, paths_df$V6, 
  paths_df$V7, paths_df$V8, paths_df$V9, paths_df$V10, paths_df$V11, paths_df$V12, paths_df$V13))
node_to_G30 = node_to_G30[!is.na(node_to_G30)]
# summarize results
filter_edges = as.data.frame(get.edgelist(g)) %>% rename("start" = "V1", "end" = "V2") %>%
  inner_join(pat_start_end_summary) %>% 
  filter(start %in% node_to_G30 & end %in% node_to_G30)
# Plot network backbone
plot_df = filter_edges
nodes_traj = as_tibble(union(names(table(plot_df$start)), names(table(plot_df$end)))) %>% 
  rowid_to_column("id") %>% rename(label = value) %>% 
  mutate(class = substr(label, start = 1, stop = 1))
per_route_traj = plot_df %>% select(start, end, n_patient) %>% rename(weight = n_patient)
edges_traj = per_route_traj %>% 
  left_join(nodes_traj, by = c("start" = "label")) %>% rename(from = id) %>% 
  left_join(nodes_traj, by = c("end" = "label")) %>% rename(to = id) %>% select(from, to, weight)
routes_tidy_traj = tbl_graph(nodes = nodes_traj, edges = edges_traj, directed = TRUE) %>% 
  activate(edges) %>% arrange(desc(weight))

# COMMAND ----------

# plot
set.seed(1234)
# pdf(paste0("/Workspace/Users/mingzhoufu@mednet.ucla.edu/outputs/cluster/", cluster_i, "_backbone.pdf") , 
#     width = 6, height = 4)
ggraph(routes_tidy_traj, layout = 'linear', circular = TRUE) + 
  geom_node_point(aes(colour = class), size = 7) +
  geom_edge_link(aes(width = weight), alpha = 0.8, arrow = arrow(length = unit(2, "mm")), end_cap = circle(1, "mm")) + 
  scale_edge_width(range = c(0.2, 2)) + geom_node_text(aes(label = label), repel = TRUE) +
  labs(edge_width = "Number of cases") + theme_graph(base_family = 'Helvetica')
# dev.off()

# COMMAND ----------

# MAGIC %md
# MAGIC ## Backbone extraction - Modularity Vitality

# COMMAND ----------

# Function to compute modularity of a given graph
compute_modularity = function(graph) {
  community = cluster_walktrap(graph, weights = E(graph)$weight)
  modularity(community)
}
# Function to compute modularity vitality of a node
modularity_vitality_node = function(g, node) {
  original_modularity = compute_modularity(g)
  g_minus_node = delete_vertices(g, node)
  new_modularity = compute_modularity(g_minus_node)
  return (original_modularity - new_modularity)
}
# Function to compute modularity vitality of an edge
modularity_vitality_edge = function(g, edge) {
  original_modularity = compute_modularity(g)
  g_minus_edge = delete_edges(g, edge)
  new_modularity = compute_modularity(g_minus_edge)
  return (original_modularity - new_modularity)
}

# COMMAND ----------

# start from significant pairs 
pat_start_end_summary_sig = pat_start_end_summary %>% 
  left_join(risk_pair_sig %>% select(exposure_ICD, outcome_ICD, coef) %>% 
  rename(start = exposure_ICD, end = outcome_ICD), by = c("start", "end")) %>% 
  filter(!is.na(coef)) %>% select(-coef) 
print(dim(pat_start_end_summary_sig)) # dim = (7796,3)
# Creating a weighted graph
g = graph_from_data_frame(data.frame(
  from = pat_start_end_summary_sig$start,
  to = pat_start_end_summary_sig$end,
  weight = pat_start_end_summary_sig$n_patient), directed = TRUE)
print(length(V(g)$name)) # 301

# COMMAND ----------

# MAGIC %md
# MAGIC ### Filter central nodes

# COMMAND ----------

# Initialize lists to store modularity vitality scores for nodes
node_vitality = numeric(vcount(g))
# Compute modularity vitality for each node
for (v in V(g)) {
  node_vitality[v] = modularity_vitality_node(g, v)
}
# create a vitality dataframe (combine nodes and edges)
node_df = data.frame(type = "node", id = V(g), vitality = abs(node_vitality)) %>% as.data.frame() %>%
  mutate(vitality = as.numeric(vitality), id = as.numeric(id))
rownames(node_df) = NULL
node_desc = data.frame(id = V(g), desc = V(g)$name) %>% as.data.frame() %>%
  mutate(id = as.numeric(id), desc = as.character(desc))
node_df = node_df %>% left_join(node_desc) %>% arrange(vitality)
save(node_df, file = paste0(raw_data_path, "clusters/full_network/node_df_0624.rda"))

# COMMAND ----------

load(file = paste0(raw_data_path, "clusters/full_network/node_df_0624.rda"))
print(dim(node_df)) # dim = (301,4)
head(node_df)

# COMMAND ----------

tail(node_df)

# COMMAND ----------

# check full modularity change with one node/edge removed at a time
g_rep = g
modularity_lst = compute_modularity(g_rep)
for (i in 1:299) {
  element_to_remove = node_df[i,]
  node_to_remove = element_to_remove$desc
  g_rep = delete_vertices(g_rep, node_to_remove)
  modularity_lst = c(modularity_lst, compute_modularity(g_rep))
}

# COMMAND ----------

mod_change_df = cbind(c(0:299), modularity_lst) %>% as.data.frame()
names(mod_change_df) = c("threshold", "modularity")
plot(mod_change_df, type = "o", col = "blue", xlab = "Step", ylab = "Modularity")

# COMMAND ----------

mod_change_df %>% filter(threshold > 250) %>% arrange(desc(modularity)) %>% head() # N = 257

# COMMAND ----------

node_df_filter = node_df[257:301,]
print(dim(node_df_filter))
node_select = node_df_filter %>% pull(desc)
print(length(node_select))

# COMMAND ----------

sort(node_select)

# COMMAND ----------

# MAGIC %md
# MAGIC ### Filter central edges

# COMMAND ----------

start_end_filter = pat_start_end_summary_sig %>% 
  filter(start %in% node_select & end %in% node_select)
print(dim(start_end_filter)) # dim = (150,3)
# Creating a weighted graph
g = graph_from_data_frame(data.frame(
  from = start_end_filter$start,
  to = start_end_filter$end,
  weight = start_end_filter$n_patient), directed = TRUE)
print(length(V(g)$name)) # 43

# COMMAND ----------

# Initialize lists to store modularity vitality scores for edges
edge_vitality = numeric(ecount(g))
# Compute modularity vitality for each edge
for (e in E(g)) {
  edge_vitality[e] = modularity_vitality_edge(g, e)
}
edge_df = data.frame(type = "edge", id = E(g), vitality = abs(edge_vitality)) %>% as.data.frame() %>%
  mutate(vitality = as.numeric(vitality), id = as.numeric(id))
edge_list = get.edgelist(g, names = TRUE)
edge_pairs = apply(edge_list, 1, function(x) paste(x[1], x[2], sep = " -> "))
edge_desc = data.frame(desc = edge_pairs, stringsAsFactors = FALSE) %>% rownames_to_column(var = "id") %>%
  mutate(id = as.numeric(id), desc = as.character(desc))
edge_df = edge_df %>% left_join(edge_desc) %>% arrange(vitality) %>%
  separate(desc, into = c("start", "end"), sep = " -> ", remove = FALSE)
print(dim(edge_df)) # dim = (150,6)
save(edge_df, file = paste0(raw_data_path, "clusters/full_network/edge_df_0624.rda"))

# COMMAND ----------

tail(edge_df)

# COMMAND ----------

# check full modularity change with one node/edge removed at a time
g_rep = g
modularity_lst = compute_modularity(g_rep)
for (i in 1:150) {
  element_to_remove = edge_df[i,]
  edge_to_remove = element_to_remove$desc
  edge_list = get.edgelist(g_rep, names = TRUE)
  edge_pairs = apply(edge_list, 1, function(x) paste(x[1], x[2], sep = " -> "))
  edge_desc = data.frame(desc = edge_pairs, stringsAsFactors = FALSE) %>% rownames_to_column(var = "id") %>%
    mutate(id = as.numeric(id), desc = as.character(desc))
  index_edge_rmv = which(edge_desc$desc == edge_to_remove)
  g_rep = delete_edges(g_rep, index_edge_rmv)
  modularity_lst = c(modularity_lst, compute_modularity(g_rep))
}

# COMMAND ----------

mod_change_df = cbind(c(0:150), modularity_lst) %>% as.data.frame()
names(mod_change_df) = c("threshold", "modularity")
plot(mod_change_df, type = "o", col = "blue", xlab = "Step", ylab = "Modularity")

# COMMAND ----------

mod_change_df %>% filter(threshold > 99 & threshold < 140) %>% arrange(desc(modularity)) %>% head(10) # N = 252

# COMMAND ----------

edge_df_filter = edge_df[100:150,]
print(dim(edge_df_filter))

# COMMAND ----------

c(edge_df_filter$start, edge_df_filter$end) %>% unique() %>% length()

# COMMAND ----------

g_rep2 = g
for (i in 1:100) {
  element_to_remove = edge_df[i,]
  edge_to_remove = element_to_remove$desc
  edge_list = get.edgelist(g_rep2, names = TRUE)
  edge_pairs = apply(edge_list, 1, function(x) paste(x[1], x[2], sep = " -> "))
  edge_desc = data.frame(desc = edge_pairs, stringsAsFactors = FALSE) %>% rownames_to_column(var = "id") %>%
    mutate(id = as.numeric(id), desc = as.character(desc))
  index_edge_rmv = which(edge_desc$desc == edge_to_remove)
  g_rep2 = delete_edges(g_rep2, index_edge_rmv)
}
compute_modularity(g_rep2)

# COMMAND ----------

# Get all paths from all vertices to vertex G30
allPathsTo_G30 = lapply(1:length(V(g_rep2)), function(v) {
  if(V(g_rep2)$name[v] != "G30") { # Ensure not to find paths from G30 to G30
    paths = all_simple_paths(g_rep2, from = V(g_rep2)$name[v], to = "G30")
    return(lapply(paths, function(p) V(g_rep2)[p]$name))
  }
})
allPathsTo_G30 = Filter(function(x) !is.null(x) && length(x) > 0, allPathsTo_G30)
allPathsTo_G30 = unlist(allPathsTo_G30, recursive = FALSE)
# Find the maximum length of the paths
max_length = max(sapply(allPathsTo_G30, length))
# Pad the shorter paths with NAs to make all paths equal in length
allPathsTo_G30 = lapply(allPathsTo_G30, function(x) {
  length(x) = max_length  
  return(x)
})
# Convert the list of paths to a data frame
paths_df = as.data.frame(do.call(rbind, allPathsTo_G30), stringsAsFactors = FALSE) %>% filter(!is.na(V3))
print(dim(paths_df)) # dim = (27260,23)
head(paths_df)

# COMMAND ----------

# plot
set.seed(1234)
# pdf(paste0("/Workspace/Users/mingzhoufu@mednet.ucla.edu/outputs/cluster/", cluster_i, "_backbone.pdf") , 
#     width = 6, height = 4)
ggraph(routes_tidy_traj, layout = 'graphopt') + 
  geom_node_point(aes(colour = class), size = 7) +
  geom_edge_link(aes(width = weight), alpha = 0.8, arrow = arrow(length = unit(2, "mm")), end_cap = circle(1, "mm")) + 
  scale_edge_width(range = c(0.2, 2)) + geom_node_text(aes(label = label), repel = TRUE) +
  labs(edge_width = "Number of cases") + theme_graph(base_family = 'Helvetica')
# dev.off()

# COMMAND ----------

# MAGIC %md
# MAGIC # 2. Clustered Networks

# COMMAND ----------

head(pat_traj_final_cluster)

# COMMAND ----------

pat_cluster = pat_traj_final_cluster %>% select(PatientID, Kmeans) %>% unique()
table(pat_cluster$Kmeans)

# COMMAND ----------

# MAGIC %md
# MAGIC ## Cluster 1: Mental disorders trajectory

# COMMAND ----------

cluster_i = 1

# COMMAND ----------

# MAGIC %md
# MAGIC ### Raw network

# COMMAND ----------

pat_cluster_i_long = pat_traj_final_cluster %>% filter(Kmeans == cluster_i) %>% 
  select(PatientID, paste0("step", c(1:9))) %>% 
  pivot_longer(cols = paste0("step", c(1:9)), names_to = "step", values_to = "ICD_3digit") %>% drop_na() %>% 
  mutate(step = as.numeric(str_sub(step, start = 5, end = -1L))) 
id_lst = unique(pat_cluster_i_long$PatientID)
print(length(id_lst)) # 1448
for (n in 1:length(id_lst)) {
  patient_id = id_lst[n]
  # print(paste0("Patinet ", n, ": ", patient_id))
  pat_rec = pat_cluster_i_long %>% filter(PatientID == patient_id) %>% 
    group_by(step) %>% mutate(traj = row_number()) %>% ungroup()
  n_traj = max(pat_rec$traj)
  for (i in 1:n_traj) {
    icd_vec = pat_rec %>% filter(traj == i) %>% pull(ICD_3digit)
    start_end_sub = data.frame(start = icd_vec[-length(icd_vec)], end = icd_vec[-1]) %>% 
      mutate(PatientID = patient_id, traj = i, n = 1)
    if (i == 1) {patient_start_end = start_end_sub}
    else {patient_start_end = rbind(patient_start_end, start_end_sub)}
  }
  if (n == 1) {start_end_full = patient_start_end}
  else {start_end_full = rbind(start_end_full, patient_start_end)}
}
pat_start_end_summary = start_end_full %>% select(-traj) %>% unique() %>% 
  group_by(start, end) %>% summarise(n_patient = n()) %>% filter(start != "G30") %>% ungroup()
print(dim(pat_start_end_summary)) # dim = (1848,3)

# COMMAND ----------

pat_cluster_i_long %>% filter(ICD_3digit == "F32") %>% pull(PatientID) %>% unique() %>% length()

# COMMAND ----------

# filter based on significant pairs and #pat followed (>=0.5%)
pat_start_end_summary_sig = pat_start_end_summary %>% 
  left_join(risk_pair_sig %>% select(exposure_ICD, outcome_ICD, coef) %>% 
  rename(start = exposure_ICD, end = outcome_ICD), by = c("start", "end")) %>% 
  filter(!is.na(coef)) %>% select(-coef) %>% filter(n_patient >= length(id_lst)*0.005)
print(dim(pat_start_end_summary_sig)) # dim = (87,3)
# Creating a weighted graph
g = graph_from_data_frame(data.frame(
  from = pat_start_end_summary_sig$start,
  to = pat_start_end_summary_sig$end,
  weight = pat_start_end_summary_sig$n_patient), directed = TRUE)
print(length(V(g)$name)) # 38
# Get all paths from all vertices to vertex G30 (remove edges that cannot lead to G30 b/c we only care about paths to G30)
allPathsTo_G30 = lapply(1:length(V(g)), function(v) {
  if(V(g)$name[v] != "G30") { # Ensure not to find paths from G30 to G30
    paths = all_simple_paths(g, from = V(g)$name[v], to = "G30")
    return(lapply(paths, function(p) V(g)[p]$name))
  }
})
allPathsTo_G30 = Filter(function(x) !is.null(x) && length(x) > 0, allPathsTo_G30)
allPathsTo_G30 = unlist(allPathsTo_G30, recursive = FALSE)
# Find the maximum length of the paths
max_length = max(sapply(allPathsTo_G30, length))
# Pad the shorter paths with NAs to make all paths equal in length
allPathsTo_G30 = lapply(allPathsTo_G30, function(x) {
  length(x) = max_length  
  return(x)
})
# Convert the list of paths to a data frame
paths_df = as.data.frame(do.call(rbind, allPathsTo_G30), stringsAsFactors = FALSE) %>% filter(!is.na(V3))
print(dim(paths_df)) # dim = (24082,13)

# COMMAND ----------

# restrict to nodes to G30
node_to_G30 = unique(c(paths_df$V1, paths_df$V2, paths_df$V3, paths_df$V4, paths_df$V5, paths_df$V6, 
  paths_df$V7, paths_df$V8, paths_df$V9, paths_df$V10, paths_df$V11, paths_df$V12, paths_df$V13))
node_to_G30 = node_to_G30[!is.na(node_to_G30)]
# summarize results
filter_edges = as.data.frame(get.edgelist(g)) %>% rename("start" = "V1", "end" = "V2") %>%
  inner_join(pat_start_end_summary) %>% 
  filter(start %in% node_to_G30 & end %in% node_to_G30)
# Plot network backbone
plot_df = filter_edges
nodes_traj = as_tibble(union(names(table(plot_df$start)), names(table(plot_df$end)))) %>% 
  rowid_to_column("id") %>% rename(label = value) %>% 
  mutate(class = substr(label, start = 1, stop = 1))
per_route_traj = plot_df %>% select(start, end, n_patient) %>% rename(weight = n_patient)
edges_traj = per_route_traj %>% 
  left_join(nodes_traj, by = c("start" = "label")) %>% rename(from = id) %>% 
  left_join(nodes_traj, by = c("end" = "label")) %>% rename(to = id) %>% select(from, to, weight)
routes_tidy_traj = tbl_graph(nodes = nodes_traj, edges = edges_traj, directed = TRUE) %>% 
  activate(edges) %>% arrange(desc(weight))

# COMMAND ----------

# plot
set.seed(1234)
# pdf(paste0("/Workspace/Users/mingzhoufu@mednet.ucla.edu/outputs/cluster/", cluster_i, "_backbone.pdf") , 
#     width = 6, height = 4)
ggraph(routes_tidy_traj, layout = 'graphopt') + 
  geom_node_point(aes(colour = class), size = 7) +
  geom_edge_link(aes(width = weight), alpha = 0.8, arrow = arrow(length = unit(2, "mm")), end_cap = circle(1, "mm")) + 
  scale_edge_width(range = c(0.2, 2)) + geom_node_text(aes(label = label), repel = TRUE) +
  labs(edge_width = "Number of cases") + theme_graph(base_family = 'Helvetica')
# dev.off()

# COMMAND ----------

# MAGIC %md
# MAGIC ### Backbone extraction - Modularity Vitality

# COMMAND ----------

# start from significant pairs 
pat_start_end_summary_sig = pat_start_end_summary %>% 
  left_join(risk_pair_sig %>% select(exposure_ICD, outcome_ICD, coef) %>% 
  rename(start = exposure_ICD, end = outcome_ICD), by = c("start", "end")) %>% 
  filter(!is.na(coef)) %>% select(-coef) %>% filter(n_patient >= length(id_lst)*0.005)
print(dim(pat_start_end_summary_sig)) # dim = (87,3)
# Creating a weighted graph
g = graph_from_data_frame(data.frame(
  from = pat_start_end_summary_sig$start,
  to = pat_start_end_summary_sig$end,
  weight = pat_start_end_summary_sig$n_patient), directed = TRUE)
print(length(V(g)$name)) # 38

# COMMAND ----------

# MAGIC %md
# MAGIC #### Sort modularity

# COMMAND ----------

# Initialize lists to store modularity vitality scores for nodes
node_vitality = numeric(vcount(g))
# Compute modularity vitality for each node
for (v in V(g)) {
  node_vitality[v] = modularity_vitality_node(g, v)
}
# create a vitality dataframe (combine nodes and edges)
node_df = data.frame(type = "node", id = V(g), vitality = abs(node_vitality)) %>% as.data.frame() %>%
  mutate(vitality = as.numeric(vitality), id = as.numeric(id))
rownames(node_df) = NULL
node_desc = data.frame(id = V(g), desc = V(g)$name) %>% as.data.frame() %>%
  mutate(id = as.numeric(id), desc = as.character(desc))
node_df = node_df %>% left_join(node_desc) %>% arrange(vitality)
print(dim(node_df)) # dim = (38,4)
tail(node_df)

# COMMAND ----------

# Initialize lists to store modularity vitality scores for edges
edge_vitality = numeric(ecount(g))
# Compute modularity vitality for each edge
for (e in E(g)) {
  edge_vitality[e] = modularity_vitality_edge(g, e)
}
edge_df = data.frame(type = "edge", id = E(g), vitality = abs(edge_vitality)) %>% as.data.frame() %>%
  mutate(vitality = as.numeric(vitality), id = as.numeric(id))
edge_list = get.edgelist(g, names = TRUE)
edge_pairs = apply(edge_list, 1, function(x) paste(x[1], x[2], sep = " -> "))
edge_desc = data.frame(desc = edge_pairs, stringsAsFactors = FALSE) %>% rownames_to_column(var = "id") %>%
  mutate(id = as.numeric(id), desc = as.character(desc))
edge_df = edge_df %>% left_join(edge_desc) %>% arrange(vitality) %>%
  separate(desc, into = c("start", "end"), sep = " -> ", remove = FALSE)
print(dim(edge_df)) # dim = (87,6)
print(tail(edge_df))

# COMMAND ----------

edge_df_short = edge_df %>% select(-c(start, end))
combined_modularity_df = rbind(node_df, edge_df_short) %>% as.data.frame() %>% arrange(vitality)
print(dim(combined_modularity_df)) # dim = (125,4)
print(tail(combined_modularity_df))

# COMMAND ----------

# check full modularity change with one node/edge removed at a time
g_rep = g
modularity_lst = compute_modularity(g_rep)
for (i in 1:123) {
  element_to_remove = combined_modularity_df[i,]
  if (element_to_remove$type == "node") {
    node_to_remove = element_to_remove$desc
    g_rep = delete_vertices(g_rep, node_to_remove)
  } else if (element_to_remove$type == "edge") {
    edge_to_remove = element_to_remove$desc
    edge_list = get.edgelist(g_rep, names = TRUE)
    edge_pairs = apply(edge_list, 1, function(x) paste(x[1], x[2], sep = " -> "))
    edge_desc = data.frame(desc = edge_pairs, stringsAsFactors = FALSE) %>% rownames_to_column(var = "id") %>%
      mutate(id = as.numeric(id), desc = as.character(desc))
    index_edge_rmv = which(edge_desc$desc == edge_to_remove)
    g_rep = delete_edges(g_rep, index_edge_rmv)
  }
  modularity_lst = c(modularity_lst, compute_modularity(g_rep))
}
mod_change_df = cbind(c(0:123), modularity_lst) %>% as.data.frame()
names(mod_change_df) = c("threshold", "modularity")
plot(mod_change_df, type = "o", col = "blue", xlab = "Step", ylab = "Modularity")

# COMMAND ----------

mod_change_df %>% filter(threshold < 100) %>% arrange(desc(modularity)) %>% head() # N = 69

# COMMAND ----------

node_edge_filter = combined_modularity_df[69:125,]
print(dim(node_edge_filter)) # dim = (57,4)
print(tail(node_edge_filter))
node_edge_select = node_edge_filter %>% pull(desc)
print(table(node_edge_filter$type))
edge_df_filter = edge_df %>% filter(desc %in% node_edge_filter$desc)
edge_filter_node = unique(c(edge_df_filter$start, edge_df_filter$end))
node_filter = node_edge_filter %>% filter(type == "node") %>% pull(desc) %>% unique()
final_node = c("G30", intersect(node_filter, edge_filter_node)) %>% unique()
print(length(final_node)) # 17
edge_df_final = edge_df_filter %>% filter(start %in% final_node & end %in% final_node) %>% left_join(pat_start_end_summary_sig) 
print(dim(edge_df_final)) # dim = (29,7)

# COMMAND ----------

# Creating a weighted graph
g_filtered = graph_from_data_frame(data.frame(
  from = edge_df_final$start,
  to = edge_df_final$end,
  weight = edge_df_final$n_patient), directed = TRUE)
print(length(V(g_filtered)$name)) # 17
# Get all paths from all vertices to vertex G30 (remove edges that cannot lead to G30 b/c we only care about paths to G30)
allPathsTo_G30 = lapply(1:length(V(g_filtered)), function(v) {
  if(V(g_filtered)$name[v] != "G30") { # Ensure not to find paths from G30 to G30
    paths = all_simple_paths(g_filtered, from = V(g_filtered)$name[v], to = "G30")
    return(lapply(paths, function(p) V(g_filtered)[p]$name))
  }
})
allPathsTo_G30 = Filter(function(x) !is.null(x) && length(x) > 0, allPathsTo_G30)
allPathsTo_G30 = unlist(allPathsTo_G30, recursive = FALSE)
# Find the maximum length of the paths
max_length = max(sapply(allPathsTo_G30, length))
# Pad the shorter paths with NAs to make all paths equal in length
allPathsTo_G30 = lapply(allPathsTo_G30, function(x) {
  length(x) = max_length  
  return(x)
})
# Convert the list of paths to a data frame
paths_df = as.data.frame(do.call(rbind, allPathsTo_G30), stringsAsFactors = FALSE) %>% filter(!is.na(V3))
print(dim(paths_df)) # dim = (164,8)

# COMMAND ----------

# restrict to nodes to G30
node_to_G30 = unique(c(paths_df$V1, paths_df$V2, paths_df$V3, paths_df$V4, paths_df$V5, paths_df$V6, paths_df$V7, paths_df$V8))
node_to_G30 = node_to_G30[!is.na(node_to_G30)]
# summarize results
filter_edges = as.data.frame(get.edgelist(g)) %>% rename("start" = "V1", "end" = "V2") %>%
  inner_join(pat_start_end_summary) %>% 
  filter(start %in% node_to_G30 & end %in% node_to_G30)

# COMMAND ----------

filter_edges %>% arrange(n_patient) %>% head(20)

# COMMAND ----------

filter_edges %>% filter(end == "F32") %>% arrange(desc(n_patient))

# COMMAND ----------

download_path = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/output/2024_09_16/"
# save(filter_edges, file = paste0(download_path, "filter_edges_", cluster_i, ".rda"))
write.csv(filter_edges, file = paste0(download_path, "filter_edges_", cluster_i, ".csv"), row.names = FALSE)

# COMMAND ----------

save(filter_edges, file = paste0(raw_data_path, "clusters/filter_edges_", cluster_i, ".rda"))

# COMMAND ----------

# Plot network backbone
plot_df = filter_edges
nodes_traj = as_tibble(union(names(table(plot_df$start)), names(table(plot_df$end)))) %>% 
  rowid_to_column("id") %>% rename(label = value) %>% 
  mutate(class = substr(label, start = 1, stop = 1))
per_route_traj = plot_df %>% select(start, end, n_patient) %>% rename(weight = n_patient)
edges_traj = per_route_traj %>% 
  left_join(nodes_traj, by = c("start" = "label")) %>% rename(from = id) %>% 
  left_join(nodes_traj, by = c("end" = "label")) %>% rename(to = id) %>% select(from, to, weight)
routes_tidy_traj = tbl_graph(nodes = nodes_traj, edges = edges_traj, directed = TRUE) %>% 
  activate(edges) %>% arrange(desc(weight))
# plot
set.seed(1347)
# pdf(paste0("/Workspace/Users/mingzhoufu@mednet.ucla.edu/outputs/cluster/", cluster_i, "_backbone.pdf") , 
#     width = 6, height = 4)
ggraph(routes_tidy_traj, layout = 'graphopt') + 
  geom_node_point(aes(colour = class), size = 14) +
  geom_edge_link(aes(width = weight), alpha = 0.8, arrow = arrow(length = unit(2, "mm")), end_cap = circle(2, "mm")) + 
  scale_edge_width(range = c(0.2, 2)) + geom_node_text(aes(label = label), repel = TRUE, size = 6) +
  labs(edge_width = "Number of cases") + theme_graph(base_family = 'Calibri')
# dev.off()

# COMMAND ----------

# MAGIC %md
# MAGIC #### Central nodes

# COMMAND ----------

central_nodes = combined_modularity_df %>% filter(desc %in% final_node & desc != "G30") %>% arrange(desc(vitality)) %>% head(5) %>% pull(desc)
# Add descriptions
icd10cm2016_3digit = icd10cm2016 %>% filter(nchar(code) == 3) %>% 
  mutate(code = as.character(code))
central_nodes_df = cbind(c(1:5), central_nodes) %>% as.data.frame() %>% 
  left_join(icd10cm2016_3digit, by = c("central_nodes" = "three_digit")) %>% 
  select(central_nodes, long_desc, sub_chapter, chapter)
central_nodes_df

# COMMAND ----------

# MAGIC %md
# MAGIC #### Trajectories in backbones

# COMMAND ----------

load(file = paste0(raw_data_path, "trajectory_clean/final/impute_wide_3step_0622.rda"))
print(dim(impute_wide_3step)) # dim = (6591,4)

# COMMAND ----------

# We will use their full trajectories (without greedy concatenation, i.e., subsets of their full trajectories)
id_lst = unique(pat_cluster_i_long$PatientID)
print(length(id_lst)) # 1448
# combine full results
folder_path = paste0(raw_data_path, "trajectory_clean/no_impute/short_traj/pat_full/") 
# List all .rda files from the folder
files = list.files(path = folder_path, pattern = "\\.rda$")
print(length(files)) # 4191

# COMMAND ----------

num_files = sapply(strsplit(files, "_"), `[`, 1) %>% as.numeric() %>% sort()
to_load_files = c()
for (id in id_lst) {
  # Use grep to find elements containing patient_id
  matching_file = files[grep(paste0("_", id, ".rda"), files)]
  to_load_files = c(to_load_files, matching_file) %>% unique()
}
print(length(to_load_files)) # should be exact same as length(id_lst)
# load in their full files
for (i in 1:length(to_load_files)){
  print(i)
  load(paste0(folder_path, to_load_files[i]))
  paths_pat_full_final$PatientID = as.character(paths_pat_full_final$PatientID)
  print(dim(paths_pat_full_final))
  if (i == 1) {
    paths_pat_combine_full = paths_pat_full_final
  } else {
    paths_pat_combine_full = bind_rows(paths_pat_combine_full, paths_pat_full_final) %>% unique()
  }
}
print(dim(paths_pat_combine_full)) # dim = (623156,17)
names(impute_wide_3step)= c("PatientID", "V1", "V2", "V3")
paths_pat_combine_cluster = bind_rows(paths_pat_combine_full, impute_wide_3step)
print(dim(paths_pat_combine_cluster)) # dim = (629747,17)
save(paths_pat_combine_cluster, file = paste0(raw_data_path, "trajectory_clean/final/paths_pat_combine_cluster_", cluster_i, ".rda"))

# COMMAND ----------

load(file = paste0(raw_data_path, "trajectory_clean/final/paths_pat_combine_cluster_", cluster_i, ".rda"))
pat_cluster_i_wide_full = paths_pat_combine_cluster %>% filter(!is.na(V3)) %>% unique() %>% 
  group_by(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16) %>% 
  summarise(n_patient = length(unique(PatientID))) %>% ungroup() 
pat_cluster_i_wide_full_short = pat_cluster_i_wide_full %>% filter(is.na(V9)) %>%
  select(V1, V2, V3, V4, V5, V6, V7, V8, n_patient)
final_pat_count = pat_cluster_i_wide_full_short %>% inner_join(paths_df) %>% arrange(desc(n_patient))
head(final_pat_count, 20)

# COMMAND ----------

# MAGIC %md
# MAGIC #### Look at >3-digit ICDs

# COMMAND ----------

load(file = paste0(raw_data_path, "patient_data/final/AD_cluster_ICD_filter_0624.rda"))
head(AD_cluster_ICD_filter)

# COMMAND ----------

length(id_lst)

# COMMAND ----------

test_icd = "F32"
pat_detail_diag = AD_cluster_ICD_filter %>% filter(PatientID %in% id_lst & ICD_3digit == test_icd) %>% 
  group_by(ICD) %>% summarise(freq = length(unique(PatientID))) 
pie_desc = pat_detail_diag %>% 
  # remove "." in the code
  mutate(code = str_remove(ICD, "\\.")) %>% left_join(icd10cm2016_short, by = c("code" = "code")) %>% 
  drop_na() %>% select(-code) %>% filter(ICD != test_icd) %>% mutate(prop = round(freq / length(id_lst)*100, 2)) %>%
  mutate(freq_prop = paste0(freq, " (", prop, "%)")) %>% select(ICD, long_desc, freq_prop)
pie_desc

# COMMAND ----------

tt = AD_cluster_ICD_filter %>% filter(PatientID %in% id_lst & ICD_3digit == test_icd) %>% pull(PatientID) %>% unique() %>% length()
pat_detail_diag$freq/tt*100

# COMMAND ----------

pie_desc = pat_detail_diag %>% 
  # remove "." in the code
  mutate(code = str_remove(ICD, "\\.")) %>% left_join(icd10cm2016_short, by = c("code" = "code")) %>% 
  drop_na() %>% select(-code) %>% filter(ICD != test_icd) %>% mutate(prop = round(freq / 1392, 2)) %>%
  mutate(freq_prop = paste0(freq, " (", prop, "%)")) %>% select(ICD, long_desc, freq_prop)
pie_desc

# COMMAND ----------

pie_desc$freq_prop

# COMMAND ----------

# MAGIC %md
# MAGIC ## Cluster 2: Cerebral vascular disease trajectory

# COMMAND ----------

cluster_i = 2

# COMMAND ----------

# MAGIC %md
# MAGIC ### Raw network

# COMMAND ----------

pat_cluster_i_long = pat_traj_final_cluster %>% filter(Kmeans == cluster_i) %>% 
  select(PatientID, paste0("step", c(1:9))) %>% 
  pivot_longer(cols = paste0("step", c(1:9)), names_to = "step", values_to = "ICD_3digit") %>% drop_na() %>% 
  mutate(step = as.numeric(str_sub(step, start = 5, end = -1L))) 
id_lst = unique(pat_cluster_i_long$PatientID)
print(length(id_lst)) # 1446
for (n in 1:length(id_lst)) {
  patient_id = id_lst[n]
  # print(paste0("Patinet ", n, ": ", patient_id))
  pat_rec = pat_cluster_i_long %>% filter(PatientID == patient_id) %>% 
    group_by(step) %>% mutate(traj = row_number()) %>% ungroup()
  n_traj = max(pat_rec$traj)
  for (i in 1:n_traj) {
    icd_vec = pat_rec %>% filter(traj == i) %>% pull(ICD_3digit)
    start_end_sub = data.frame(start = icd_vec[-length(icd_vec)], end = icd_vec[-1]) %>% 
      mutate(PatientID = patient_id, traj = i, n = 1)
    if (i == 1) {patient_start_end = start_end_sub}
    else {patient_start_end = rbind(patient_start_end, start_end_sub)}
  }
  if (n == 1) {start_end_full = patient_start_end}
  else {start_end_full = rbind(start_end_full, patient_start_end)}
}
pat_start_end_summary = start_end_full %>% select(-traj) %>% unique() %>% 
  group_by(start, end) %>% summarise(n_patient = n()) %>% filter(start != "G30") %>% ungroup()
print(dim(pat_start_end_summary)) # dim = (6619,3)

# COMMAND ----------

pat_cluster_i_long %>% filter(ICD_3digit == "I67") %>% pull(PatientID) %>% unique() %>% length()

# COMMAND ----------

# filter based on significant pairs and #pat followed (>=0.5%)
pat_start_end_summary_sig = pat_start_end_summary %>% 
  left_join(risk_pair_sig %>% select(exposure_ICD, outcome_ICD, coef) %>% 
  rename(start = exposure_ICD, end = outcome_ICD), by = c("start", "end")) %>% 
  filter(!is.na(coef)) %>% select(-coef) %>% filter(n_patient >= length(id_lst)*0.005)
print(dim(pat_start_end_summary_sig)) # dim = (202,3)
# Creating a weighted graph
g = graph_from_data_frame(data.frame(
  from = pat_start_end_summary_sig$start,
  to = pat_start_end_summary_sig$end,
  weight = pat_start_end_summary_sig$n_patient), directed = TRUE)
print(length(V(g)$name)) # 79
# Get all paths from all vertices to vertex G30 (remove edges that cannot lead to G30 b/c we only care about paths to G30)
allPathsTo_G30 = lapply(1:length(V(g)), function(v) {
  if(V(g)$name[v] != "G30") { # Ensure not to find paths from G30 to G30
    paths = all_simple_paths(g, from = V(g)$name[v], to = "G30")
    return(lapply(paths, function(p) V(g)[p]$name))
  }
})
allPathsTo_G30 = Filter(function(x) !is.null(x) && length(x) > 0, allPathsTo_G30)
allPathsTo_G30 = unlist(allPathsTo_G30, recursive = FALSE)
# Find the maximum length of the paths
max_length = max(sapply(allPathsTo_G30, length))
# Pad the shorter paths with NAs to make all paths equal in length
allPathsTo_G30 = lapply(allPathsTo_G30, function(x) {
  length(x) = max_length  
  return(x)
})
# Convert the list of paths to a data frame
paths_df = as.data.frame(do.call(rbind, allPathsTo_G30), stringsAsFactors = FALSE) %>% filter(!is.na(V3))
print(dim(paths_df)) # dim = (9454203,24)

# COMMAND ----------

# restrict to nodes to G30
node_to_G30 = unique(c(paths_df$V1, paths_df$V2, paths_df$V3, paths_df$V4, paths_df$V5, paths_df$V6, 
  paths_df$V7, paths_df$V8, paths_df$V9, paths_df$V10, paths_df$V11, paths_df$V12, paths_df$V13, paths_df$V14, 
  paths_df$V15, paths_df$V16, paths_df$V17, paths_df$V18, paths_df$V19, paths_df$V20, paths_df$V21, paths_df$V22, 
  paths_df$V23, paths_df$V24))
node_to_G30 = node_to_G30[!is.na(node_to_G30)]
# summarize results
filter_edges = as.data.frame(get.edgelist(g)) %>% rename("start" = "V1", "end" = "V2") %>%
  inner_join(pat_start_end_summary) %>% 
  filter(start %in% node_to_G30 & end %in% node_to_G30)
# Plot network backbone
plot_df = filter_edges
nodes_traj = as_tibble(union(names(table(plot_df$start)), names(table(plot_df$end)))) %>% 
  rowid_to_column("id") %>% rename(label = value) %>% 
  mutate(class = substr(label, start = 1, stop = 1))
per_route_traj = plot_df %>% select(start, end, n_patient) %>% rename(weight = n_patient)
edges_traj = per_route_traj %>% 
  left_join(nodes_traj, by = c("start" = "label")) %>% rename(from = id) %>% 
  left_join(nodes_traj, by = c("end" = "label")) %>% rename(to = id) %>% select(from, to, weight)
routes_tidy_traj = tbl_graph(nodes = nodes_traj, edges = edges_traj, directed = TRUE) %>% 
  activate(edges) %>% arrange(desc(weight))

# COMMAND ----------

# plot
set.seed(1234)
# pdf(paste0("/Workspace/Users/mingzhoufu@mednet.ucla.edu/outputs/cluster/", cluster_i, "_backbone.pdf") , 
#     width = 6, height = 4)
ggraph(routes_tidy_traj, layout = 'graphopt') + 
  geom_node_point(aes(colour = class), size = 7) +
  geom_edge_link(aes(width = weight), alpha = 0.8, arrow = arrow(length = unit(2, "mm")), end_cap = circle(1, "mm")) + 
  scale_edge_width(range = c(0.2, 2)) + geom_node_text(aes(label = label), repel = TRUE) +
  labs(edge_width = "Number of cases") + theme_graph(base_family = 'Helvetica')
# dev.off()

# COMMAND ----------

# MAGIC %md
# MAGIC ### Backbone extraction - Modularity Vitality

# COMMAND ----------

# start from significant pairs 
pat_start_end_summary_sig = pat_start_end_summary %>% 
  left_join(risk_pair_sig %>% select(exposure_ICD, outcome_ICD, coef) %>% 
  rename(start = exposure_ICD, end = outcome_ICD), by = c("start", "end")) %>% 
  filter(!is.na(coef)) %>% select(-coef) %>% filter(n_patient >= length(id_lst)*0.005)
print(dim(pat_start_end_summary_sig)) # dim = (202,3)
# Creating a weighted graph
g = graph_from_data_frame(data.frame(
  from = pat_start_end_summary_sig$start,
  to = pat_start_end_summary_sig$end,
  weight = pat_start_end_summary_sig$n_patient), directed = TRUE)
print(length(V(g)$name)) # 79

# COMMAND ----------

# MAGIC %md
# MAGIC #### Sort modularity

# COMMAND ----------

# Initialize lists to store modularity vitality scores for nodes
node_vitality = numeric(vcount(g))
# Compute modularity vitality for each node
for (v in V(g)) {
  node_vitality[v] = modularity_vitality_node(g, v)
}
# create a vitality dataframe (combine nodes and edges)
node_df = data.frame(type = "node", id = V(g), vitality = abs(node_vitality)) %>% as.data.frame() %>%
  mutate(vitality = as.numeric(vitality), id = as.numeric(id))
rownames(node_df) = NULL
node_desc = data.frame(id = V(g), desc = V(g)$name) %>% as.data.frame() %>%
  mutate(id = as.numeric(id), desc = as.character(desc))
node_df = node_df %>% left_join(node_desc) %>% arrange(vitality)
print(dim(node_df)) # dim = (79,4)
tail(node_df)

# COMMAND ----------

# Initialize lists to store modularity vitality scores for edges
edge_vitality = numeric(ecount(g))
# Compute modularity vitality for each edge
for (e in E(g)) {
  edge_vitality[e] = modularity_vitality_edge(g, e)
}
edge_df = data.frame(type = "edge", id = E(g), vitality = abs(edge_vitality)) %>% as.data.frame() %>%
  mutate(vitality = as.numeric(vitality), id = as.numeric(id))
edge_list = get.edgelist(g, names = TRUE)
edge_pairs = apply(edge_list, 1, function(x) paste(x[1], x[2], sep = " -> "))
edge_desc = data.frame(desc = edge_pairs, stringsAsFactors = FALSE) %>% rownames_to_column(var = "id") %>%
  mutate(id = as.numeric(id), desc = as.character(desc))
edge_df = edge_df %>% left_join(edge_desc) %>% arrange(vitality) %>%
  separate(desc, into = c("start", "end"), sep = " -> ", remove = FALSE)
print(dim(edge_df)) # dim = (202,6)
print(tail(edge_df))

# COMMAND ----------

edge_df_short = edge_df %>% select(-c(start, end))
combined_modularity_df = rbind(node_df, edge_df_short) %>% as.data.frame() %>% arrange(vitality)
print(dim(combined_modularity_df)) # dim = (281,4)
print(tail(combined_modularity_df))

# COMMAND ----------

which(combined_modularity_df$desc == "G30")

# COMMAND ----------

# check full modularity change with one node/edge removed at a time
g_rep = g
modularity_lst = compute_modularity(g_rep)
for (i in 1:280) {
  element_to_remove = combined_modularity_df[i,]
  if (element_to_remove$type == "node") {
    node_to_remove = element_to_remove$desc
    g_rep = delete_vertices(g_rep, node_to_remove)
  } else if (element_to_remove$type == "edge") {
    edge_to_remove = element_to_remove$desc
    edge_list = get.edgelist(g_rep, names = TRUE)
    edge_pairs = apply(edge_list, 1, function(x) paste(x[1], x[2], sep = " -> "))
    edge_desc = data.frame(desc = edge_pairs, stringsAsFactors = FALSE) %>% rownames_to_column(var = "id") %>%
      mutate(id = as.numeric(id), desc = as.character(desc))
    index_edge_rmv = which(edge_desc$desc == edge_to_remove)
    g_rep = delete_edges(g_rep, index_edge_rmv)
  }
  modularity_lst = c(modularity_lst, compute_modularity(g_rep))
}
mod_change_df = cbind(c(0:280), modularity_lst) %>% as.data.frame()
names(mod_change_df) = c("threshold", "modularity")
plot(mod_change_df, type = "o", col = "blue", xlab = "Step", ylab = "Modularity")

# COMMAND ----------

mod_change_df %>% filter(threshold > 200) %>% arrange(desc(modularity)) %>% head() # N = 211

# COMMAND ----------

node_edge_filter = combined_modularity_df[211:281,]
print(dim(node_edge_filter)) # dim = (71,4)
print(tail(node_edge_filter))
node_edge_select = node_edge_filter %>% pull(desc)
print(table(node_edge_filter$type))
edge_df_filter = edge_df %>% filter(desc %in% node_edge_filter$desc)
edge_filter_node = unique(c(edge_df_filter$start, edge_df_filter$end))
node_filter = node_edge_filter %>% filter(type == "node") %>% pull(desc) %>% unique()
final_node = c("G30", intersect(node_filter, edge_filter_node)) %>% unique()
print(length(final_node)) # 25
edge_df_final = edge_df_filter %>% filter(start %in% final_node & end %in% final_node) %>% left_join(pat_start_end_summary_sig) 
print(dim(edge_df_final)) # dim = (34,7)

# COMMAND ----------

# Creating a weighted graph
g_filtered = graph_from_data_frame(data.frame(
  from = edge_df_final$start,
  to = edge_df_final$end,
  weight = edge_df_final$n_patient), directed = TRUE)
print(length(V(g_filtered)$name)) # 25
# Get all paths from all vertices to vertex G30 (remove edges that cannot lead to G30 b/c we only care about paths to G30)
allPathsTo_G30 = lapply(1:length(V(g_filtered)), function(v) {
  if(V(g_filtered)$name[v] != "G30") { # Ensure not to find paths from G30 to G30
    paths = all_simple_paths(g_filtered, from = V(g_filtered)$name[v], to = "G30")
    return(lapply(paths, function(p) V(g_filtered)[p]$name))
  }
})
allPathsTo_G30 = Filter(function(x) !is.null(x) && length(x) > 0, allPathsTo_G30)
allPathsTo_G30 = unlist(allPathsTo_G30, recursive = FALSE)
# Find the maximum length of the paths
max_length = max(sapply(allPathsTo_G30, length))
# Pad the shorter paths with NAs to make all paths equal in length
allPathsTo_G30 = lapply(allPathsTo_G30, function(x) {
  length(x) = max_length  
  return(x)
})
# Convert the list of paths to a data frame
paths_df = as.data.frame(do.call(rbind, allPathsTo_G30), stringsAsFactors = FALSE) %>% filter(!is.na(V3))
print(dim(paths_df)) # dim = (40,5)

# COMMAND ----------

# restrict to nodes to G30
node_to_G30 = unique(c(paths_df$V1, paths_df$V2, paths_df$V3, paths_df$V4, paths_df$V5))
node_to_G30 = node_to_G30[!is.na(node_to_G30)]
# summarize results
filter_edges = as.data.frame(get.edgelist(g)) %>% rename("start" = "V1", "end" = "V2") %>%
  inner_join(pat_start_end_summary) %>% 
  filter(start %in% node_to_G30 & end %in% node_to_G30)

# COMMAND ----------

filter_edges %>% arrange(n_patient) %>% head(20)

# COMMAND ----------

filter_edges %>% filter(end == "I67") %>% arrange(desc(n_patient)) %>% head()

# COMMAND ----------

download_path = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/output/2024_09_16/"
# save(filter_edges, file = paste0(download_path, "filter_edges_", cluster_i, ".rda"))
write.csv(filter_edges, file = paste0(download_path, "filter_edges_", cluster_i, ".csv"), row.names = FALSE)

# COMMAND ----------

save(filter_edges, file = paste0(raw_data_path, "clusters/filter_edges_", cluster_i, ".rda"))

# COMMAND ----------

# Plot network backbone
plot_df = filter_edges
nodes_traj = as_tibble(union(names(table(plot_df$start)), names(table(plot_df$end)))) %>% 
  rowid_to_column("id") %>% rename(label = value) %>% 
  mutate(class = substr(label, start = 1, stop = 1))
per_route_traj = plot_df %>% select(start, end, n_patient) %>% rename(weight = n_patient)
edges_traj = per_route_traj %>% 
  left_join(nodes_traj, by = c("start" = "label")) %>% rename(from = id) %>% 
  left_join(nodes_traj, by = c("end" = "label")) %>% rename(to = id) %>% select(from, to, weight)
routes_tidy_traj = tbl_graph(nodes = nodes_traj, edges = edges_traj, directed = TRUE) %>% 
  activate(edges) %>% arrange(desc(weight))
# plot
set.seed(1347)
# # pdf(paste0("/Workspace/Users/mingzhoufu@mednet.ucla.edu/outputs/cluster/", cluster_i, "_backbone.pdf") , 
# #     width = 6, height = 4)
# ggraph(routes_tidy_traj, layout = 'graphopt') + 
#   geom_node_point(aes(colour = class), size = 7) +
#   geom_edge_link(aes(width = weight), alpha = 0.8, arrow = arrow(length = unit(2, "mm")), end_cap = circle(1, "mm")) + 
#   scale_edge_width(range = c(0.2, 2)) + geom_node_text(aes(label = label), repel = TRUE) +
#   labs(edge_width = "Number of cases") + theme_graph(base_family = 'Helvetica')
# # dev.off()
ggraph(routes_tidy_traj, layout = 'graphopt') + 
  geom_node_point(aes(colour = class), size = 14) +
  geom_edge_link(aes(width = weight), alpha = 0.8, arrow = arrow(length = unit(2, "mm")), end_cap = circle(2, "mm")) + 
  scale_edge_width(range = c(0.2, 2)) + geom_node_text(aes(label = label), repel = TRUE, size = 6) +
  labs(edge_width = "Number of cases") + theme_graph(base_family = 'Calibri')

# COMMAND ----------

# MAGIC %md
# MAGIC #### Central nodes

# COMMAND ----------

central_nodes = combined_modularity_df %>% filter(desc %in% final_node & desc != "G30") %>% arrange(desc(vitality)) %>% head(5) %>% pull(desc)
# Add descriptions
icd10cm2016_3digit = icd10cm2016 %>% filter(nchar(code) == 3) %>% 
  mutate(code = as.character(code))
central_nodes_df = cbind(c(1:5), central_nodes) %>% as.data.frame() %>% 
  left_join(icd10cm2016_3digit, by = c("central_nodes" = "three_digit")) %>% 
  select(central_nodes, long_desc, sub_chapter, chapter)
central_nodes_df

# COMMAND ----------

# MAGIC %md
# MAGIC #### Trajectories in backbones

# COMMAND ----------

# We will use their full trajectories (without greedy concatenation, i.e., subsets of their full trajectories)
id_lst = unique(pat_cluster_i_long$PatientID)
print(length(id_lst)) # 1446
# combine full results
folder_path = paste0(raw_data_path, "trajectory_clean/no_impute/short_traj/pat_full/") 
# List all .rda files from the folder
files = list.files(path = folder_path, pattern = "\\.rda$")
print(length(files)) # 4191

# COMMAND ----------

num_files = sapply(strsplit(files, "_"), `[`, 1) %>% as.numeric() %>% sort()
to_load_files = c()
for (id in id_lst) {
  # Use grep to find elements containing patient_id
  matching_file = files[grep(paste0("_", id, ".rda"), files)]
  to_load_files = c(to_load_files, matching_file) %>% unique()
}
print(length(to_load_files)) # should be exact same as length(id_lst)
# load in their full files
for (i in 1:length(to_load_files)){
  print(i)
  load(paste0(folder_path, to_load_files[i]))
  paths_pat_full_final$PatientID = as.character(paths_pat_full_final$PatientID)
  print(dim(paths_pat_full_final))
  if (i == 1) {
    paths_pat_combine_full = paths_pat_full_final
  } else {
    paths_pat_combine_full = bind_rows(paths_pat_combine_full, paths_pat_full_final) %>% unique()
  }
}
print(dim(paths_pat_combine_full)) # dim = (2367334,18)

# COMMAND ----------

names(impute_wide_3step)= c("PatientID", "V1", "V2", "V3")
paths_pat_combine_cluster = bind_rows(paths_pat_combine_full, impute_wide_3step)
print(dim(paths_pat_combine_cluster)) # dim = (2373925,18)
save(paths_pat_combine_cluster, file = paste0(raw_data_path, "trajectory_clean/final/paths_pat_combine_cluster_", cluster_i, ".rda"))

# COMMAND ----------

load(file = paste0(raw_data_path, "trajectory_clean/final/paths_pat_combine_cluster_", cluster_i, ".rda"))
pat_cluster_i_wide_full = paths_pat_combine_cluster %>% filter(!is.na(V3)) %>% unique() %>% 
  group_by(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16) %>% 
  summarise(n_patient = length(unique(PatientID))) %>% ungroup() 
pat_cluster_i_wide_full_short = pat_cluster_i_wide_full %>% filter(is.na(V6)) %>%
  select(V1, V2, V3, V4, V5, n_patient)
final_pat_count = pat_cluster_i_wide_full_short %>% inner_join(paths_df) %>% arrange(desc(n_patient))
head(final_pat_count, 20)

# COMMAND ----------

# MAGIC %md
# MAGIC #### Look at >3-digit ICDs

# COMMAND ----------

length(id_lst)

# COMMAND ----------

test_icd = "I67"
pat_detail_diag = AD_cluster_ICD_filter %>% filter(PatientID %in% id_lst & ICD_3digit == test_icd) %>% 
  group_by(ICD) %>% summarise(freq = length(unique(PatientID))) 
pie_desc = pat_detail_diag %>% 
  # remove "." in the code
  mutate(code = str_remove(ICD, "\\.")) %>% left_join(icd10cm2016_short, by = c("code" = "code")) %>% 
  drop_na() %>% select(-code) %>% filter(ICD != test_icd) %>% mutate(prop = round(freq / length(id_lst)*100, 2)) %>%
  mutate(freq_prop = paste0(freq, " (", prop, "%)")) %>% select(ICD, long_desc, freq_prop)
pie_desc

# COMMAND ----------

tt = AD_cluster_ICD_filter %>% filter(PatientID %in% id_lst & ICD_3digit == test_icd) %>% pull(PatientID) %>% unique() %>% length()
pat_detail_diag$freq/tt*100

# COMMAND ----------

pat_detail_diag

# COMMAND ----------



# COMMAND ----------

# We will use their full trajectories (without greedy concatenation, i.e., subsets of their full trajectories)
folder_path = paste0(raw_data_path, "trajectory_clean/pat_full/")
id_lst = unique(pat_cluster_i_long$PatientID)
print(length(id_lst))
# List all .rda files from the folder
files = list.files(path = folder_path, pattern = "\\.rda$")
print(length(files)) # 7627
to_load_files = c()
for (id in id_lst) {
  # Use grep to find elements containing patient_id
  matching_file = files[grep(paste0("_", id, ".rda"), files)]
  to_load_files = c(to_load_files, matching_file) %>% unique()
}
print(length(to_load_files)) # should be exact same as length(id_lst)

# COMMAND ----------

# load in their full files
for (i in 1:length(to_load_files)){
  # print(i)
  load(paste0(folder_path, to_load_files[i]))
  paths_pat_full_final$PatientID = as.character(paths_pat_full_final$PatientID)
  # print(dim(paths_pat_full_final))
  if (i == 1) {
    paths_pat_combine_full = paths_pat_full_final
  } else {
    paths_pat_combine_full = bind_rows(paths_pat_combine_full, paths_pat_full_final) %>% unique()
  }
}
print(dim(paths_pat_combine_full)) # dim = (59511,18)
save(paths_pat_combine_full, file = paste0(raw_data_path, "trajectory_clean/final/paths_pat_combine_full_", cluster_i, ".rda"))

# COMMAND ----------

# load(file = paste0(raw_data_path, "trajectory_clean/final/paths_pat_combine_full_", cluster_i, ".rda"))
pat_cluster_i_wide = paths_pat_combine_full %>% filter(!is.na(V3) & is.na(V11)) %>% unique() %>% 
  group_by(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17) %>% 
  summarise(n_patient = length(unique(PatientID))) %>% ungroup() 
node_to_G30

# COMMAND ----------

# only keep trajectories include node_to_G30
traj_cluster_i = apply(pat_cluster_i_wide[, !names(pat_cluster_i_wide) %in% "n_patient"], 1, function(x) all(na.omit(x) %in% node_to_G30))
# Subset the dataframe
pat_traj_cluster_i = pat_cluster_i_wide[traj_cluster_i, ] 
# Identify columns where all values are NA
na_columns = colSums(is.na(pat_traj_cluster_i)) == nrow(pat_traj_cluster_i)
# Drop these columns
pat_traj_cluster_i = pat_traj_cluster_i[, !na_columns] %>% arrange(desc(n_patient))
pat_traj_top5 = pat_traj_cluster_i %>% head(5) 
pat_traj_top5

# COMMAND ----------

names(pat_traj_top5) = c(paste0("step", c(1:9)), "n_patient")
pat_traj_top5_i = pat_cluster_full %>% filter(cluster_pam == cluster_i) %>% 
  filter(is.na(step4)) %>% select(paste0("step", c(1:3))) %>% unique() %>% inner_join(pat_traj_top5)
pat_traj_top5_i %>% select(step1, step2, step3, n_patient) %>% arrange(desc(n_patient))

# COMMAND ----------

unique(c(pat_traj_top5$step1, pat_traj_top5$step2, pat_traj_top5$step3)) %>% as.data.frame() %>% rename("code" = ".") %>%
  left_join(icd10cm2016_3digit, by = c("code" = "three_digit")) %>% select(code, long_desc) 

# COMMAND ----------



# COMMAND ----------

test_icd = "G31"
pat_detail_diag = AD_encounters_mod %>% filter(PatientID %in% id_lst & ICD_3digit == test_icd) %>% 
  group_by(DiagnosisCode) %>% summarise(freq = length(unique(PatientID))) 
pie_desc = pat_detail_diag %>% 
  # remove "." in the code
  mutate(code = str_remove(DiagnosisCode, "\\.")) %>% left_join(icd10cm2016_short, by = c("code" = "code")) %>% 
  drop_na() %>% select(-code) %>% filter(DiagnosisCode != test_icd) %>% mutate(prop = round(freq / length(id_lst)*100, 2)) %>%
  mutate(freq_prop = paste0(freq, " (", prop, "%)")) %>% select(DiagnosisCode, long_desc, freq_prop)
pie_desc

# COMMAND ----------

pie_desc %>% select(DiagnosisCode, freq_prop)

# COMMAND ----------



# COMMAND ----------

# MAGIC %md
# MAGIC ## Cluster 3: Brain disease trajectory

# COMMAND ----------

cluster_i = 3

# COMMAND ----------

# MAGIC %md
# MAGIC #### Raw network

# COMMAND ----------

pat_cluster_i_long = pat_traj_final_cluster %>% filter(Kmeans == cluster_i) %>% 
  select(PatientID, paste0("step", c(1:9))) %>% 
  pivot_longer(cols = paste0("step", c(1:9)), names_to = "step", values_to = "ICD_3digit") %>% drop_na() %>% 
  mutate(step = as.numeric(str_sub(step, start = 5, end = -1L))) 
id_lst = unique(pat_cluster_i_long$PatientID)
print(length(id_lst)) # 3223
for (n in 1:length(id_lst)) {
  patient_id = id_lst[n]
  # print(paste0("Patinet ", n, ": ", patient_id))
  pat_rec = pat_cluster_i_long %>% filter(PatientID == patient_id) %>% 
    group_by(step) %>% mutate(traj = row_number()) %>% ungroup()
  n_traj = max(pat_rec$traj)
  for (i in 1:n_traj) {
    icd_vec = pat_rec %>% filter(traj == i) %>% pull(ICD_3digit)
    start_end_sub = data.frame(start = icd_vec[-length(icd_vec)], end = icd_vec[-1]) %>% 
      mutate(PatientID = patient_id, traj = i, n = 1)
    if (i == 1) {patient_start_end = start_end_sub}
    else {patient_start_end = rbind(patient_start_end, start_end_sub)}
  }
  if (n == 1) {start_end_full = patient_start_end}
  else {start_end_full = rbind(start_end_full, patient_start_end)}
}
pat_start_end_summary = start_end_full %>% select(-traj) %>% unique() %>% 
  group_by(start, end) %>% summarise(n_patient = n()) %>% filter(start != "G30") %>% ungroup()
print(dim(pat_start_end_summary)) # dim = (618,3)

# COMMAND ----------

pat_cluster_i_long %>% filter(ICD_3digit == "G93") %>% pull(PatientID) %>% unique() %>% length()

# COMMAND ----------

# filter based on significant pairs and #pat followed (>=0.5%)
pat_start_end_summary_sig = pat_start_end_summary %>% 
  left_join(risk_pair_sig %>% select(exposure_ICD, outcome_ICD, coef) %>% 
  rename(start = exposure_ICD, end = outcome_ICD), by = c("start", "end")) %>% 
  filter(!is.na(coef)) %>% select(-coef) %>% filter(n_patient >= length(id_lst)*0.005)
print(dim(pat_start_end_summary_sig)) # dim = (75,3)
# Creating a weighted graph
g = graph_from_data_frame(data.frame(
  from = pat_start_end_summary_sig$start,
  to = pat_start_end_summary_sig$end,
  weight = pat_start_end_summary_sig$n_patient), directed = TRUE)
print(length(V(g)$name)) # 30
# Get all paths from all vertices to vertex G30 (remove edges that cannot lead to G30 b/c we only care about paths to G30)
allPathsTo_G30 = lapply(1:length(V(g)), function(v) {
  if(V(g)$name[v] != "G30") { # Ensure not to find paths from G30 to G30
    paths = all_simple_paths(g, from = V(g)$name[v], to = "G30")
    return(lapply(paths, function(p) V(g)[p]$name))
  }
})
allPathsTo_G30 = Filter(function(x) !is.null(x) && length(x) > 0, allPathsTo_G30)
allPathsTo_G30 = unlist(allPathsTo_G30, recursive = FALSE)
# Find the maximum length of the paths
max_length = max(sapply(allPathsTo_G30, length))
# Pad the shorter paths with NAs to make all paths equal in length
allPathsTo_G30 = lapply(allPathsTo_G30, function(x) {
  length(x) = max_length  
  return(x)
})
# Convert the list of paths to a data frame
paths_df = as.data.frame(do.call(rbind, allPathsTo_G30), stringsAsFactors = FALSE) %>% filter(!is.na(V3))
print(dim(paths_df)) # dim = (1336,10)

# COMMAND ----------

# restrict to nodes to G30
node_to_G30 = unique(c(paths_df$V1, paths_df$V2, paths_df$V3, paths_df$V4, paths_df$V5, paths_df$V6, 
  paths_df$V7, paths_df$V8, paths_df$V9, paths_df$V10))
node_to_G30 = node_to_G30[!is.na(node_to_G30)]
# summarize results
filter_edges = as.data.frame(get.edgelist(g)) %>% rename("start" = "V1", "end" = "V2") %>%
  inner_join(pat_start_end_summary) %>% 
  filter(start %in% node_to_G30 & end %in% node_to_G30)
# Plot network backbone
plot_df = filter_edges
nodes_traj = as_tibble(union(names(table(plot_df$start)), names(table(plot_df$end)))) %>% 
  rowid_to_column("id") %>% rename(label = value) %>% 
  mutate(class = substr(label, start = 1, stop = 1))
per_route_traj = plot_df %>% select(start, end, n_patient) %>% rename(weight = n_patient)
edges_traj = per_route_traj %>% 
  left_join(nodes_traj, by = c("start" = "label")) %>% rename(from = id) %>% 
  left_join(nodes_traj, by = c("end" = "label")) %>% rename(to = id) %>% select(from, to, weight)
routes_tidy_traj = tbl_graph(nodes = nodes_traj, edges = edges_traj, directed = TRUE) %>% 
  activate(edges) %>% arrange(desc(weight))

# COMMAND ----------

# plot
set.seed(1234)
# pdf(paste0("/Workspace/Users/mingzhoufu@mednet.ucla.edu/outputs/cluster/", cluster_i, "_backbone.pdf") , 
#     width = 6, height = 4)
ggraph(routes_tidy_traj, layout = 'graphopt') + 
  geom_node_point(aes(colour = class), size = 7) +
  geom_edge_link(aes(width = weight), alpha = 0.8, arrow = arrow(length = unit(2, "mm")), end_cap = circle(1, "mm")) + 
  scale_edge_width(range = c(0.2, 2)) + geom_node_text(aes(label = label), repel = TRUE) +
  labs(edge_width = "Number of cases") + theme_graph(base_family = 'Helvetica')
# dev.off()

# COMMAND ----------

# MAGIC %md
# MAGIC ### Backbone extraction - Modularity Vitality

# COMMAND ----------

# start from significant pairs 
pat_start_end_summary_sig = pat_start_end_summary %>% 
  left_join(risk_pair_sig %>% select(exposure_ICD, outcome_ICD, coef) %>% 
  rename(start = exposure_ICD, end = outcome_ICD), by = c("start", "end")) %>% 
  filter(!is.na(coef)) %>% select(-coef) %>% filter(n_patient >= length(id_lst)*0.005)
print(dim(pat_start_end_summary_sig)) # dim = (75,3)
# Creating a weighted graph
g = graph_from_data_frame(data.frame(
  from = pat_start_end_summary_sig$start,
  to = pat_start_end_summary_sig$end,
  weight = pat_start_end_summary_sig$n_patient), directed = TRUE)
print(length(V(g)$name)) # 30

# COMMAND ----------

# MAGIC %md
# MAGIC #### Sort modularity

# COMMAND ----------

# Initialize lists to store modularity vitality scores for nodes
node_vitality = numeric(vcount(g))
# Compute modularity vitality for each node
for (v in V(g)) {
  node_vitality[v] = modularity_vitality_node(g, v)
}
# create a vitality dataframe (combine nodes and edges)
node_df = data.frame(type = "node", id = V(g), vitality = abs(node_vitality)) %>% as.data.frame() %>%
  mutate(vitality = as.numeric(vitality), id = as.numeric(id))
rownames(node_df) = NULL
node_desc = data.frame(id = V(g), desc = V(g)$name) %>% as.data.frame() %>%
  mutate(id = as.numeric(id), desc = as.character(desc))
node_df = node_df %>% left_join(node_desc) %>% arrange(vitality)
print(dim(node_df)) # dim = (30,4)
tail(node_df)

# COMMAND ----------

# Initialize lists to store modularity vitality scores for edges
edge_vitality = numeric(ecount(g))
# Compute modularity vitality for each edge
for (e in E(g)) {
  edge_vitality[e] = modularity_vitality_edge(g, e)
}
edge_df = data.frame(type = "edge", id = E(g), vitality = abs(edge_vitality)) %>% as.data.frame() %>%
  mutate(vitality = as.numeric(vitality), id = as.numeric(id))
edge_list = get.edgelist(g, names = TRUE)
edge_pairs = apply(edge_list, 1, function(x) paste(x[1], x[2], sep = " -> "))
edge_desc = data.frame(desc = edge_pairs, stringsAsFactors = FALSE) %>% rownames_to_column(var = "id") %>%
  mutate(id = as.numeric(id), desc = as.character(desc))
edge_df = edge_df %>% left_join(edge_desc) %>% arrange(vitality) %>%
  separate(desc, into = c("start", "end"), sep = " -> ", remove = FALSE)
print(dim(edge_df)) # dim = (75,6)
print(tail(edge_df))

# COMMAND ----------

edge_df_short = edge_df %>% select(-c(start, end))
combined_modularity_df = rbind(node_df, edge_df_short) %>% as.data.frame() %>% arrange(vitality)
print(dim(combined_modularity_df)) # dim = (105,4)
print(tail(combined_modularity_df))

# COMMAND ----------

# check full modularity change with one node/edge removed at a time
g_rep = g
modularity_lst = compute_modularity(g_rep)
for (i in 1:104) {
  element_to_remove = combined_modularity_df[i,]
  if (element_to_remove$type == "node") {
    node_to_remove = element_to_remove$desc
    g_rep = delete_vertices(g_rep, node_to_remove)
  } else if (element_to_remove$type == "edge") {
    edge_to_remove = element_to_remove$desc
    edge_list = get.edgelist(g_rep, names = TRUE)
    edge_pairs = apply(edge_list, 1, function(x) paste(x[1], x[2], sep = " -> "))
    edge_desc = data.frame(desc = edge_pairs, stringsAsFactors = FALSE) %>% rownames_to_column(var = "id") %>%
      mutate(id = as.numeric(id), desc = as.character(desc))
    index_edge_rmv = which(edge_desc$desc == edge_to_remove)
    g_rep = delete_edges(g_rep, index_edge_rmv)
  }
  modularity_lst = c(modularity_lst, compute_modularity(g_rep))
}
mod_change_df = cbind(c(0:104), modularity_lst) %>% as.data.frame()
names(mod_change_df) = c("threshold", "modularity")
plot(mod_change_df, type = "o", col = "blue", xlab = "Step", ylab = "Modularity")

# COMMAND ----------

mod_change_df %>% filter(threshold > 40) %>% arrange(desc(modularity)) %>% head() # N = 54

# COMMAND ----------

node_edge_filter = combined_modularity_df[54:105,]
print(dim(node_edge_filter)) # dim = (52,4)
print(tail(node_edge_filter))
node_edge_select = node_edge_filter %>% pull(desc)
print(table(node_edge_filter$type))
edge_df_filter = edge_df %>% filter(desc %in% node_edge_filter$desc)
edge_filter_node = unique(c(edge_df_filter$start, edge_df_filter$end))
node_filter = node_edge_filter %>% filter(type == "node") %>% pull(desc) %>% unique()
final_node = c("G30", intersect(node_filter, edge_filter_node)) %>% unique()
print(length(final_node)) # 17
edge_df_final = edge_df_filter %>% filter(start %in% final_node & end %in% final_node) %>% left_join(pat_start_end_summary_sig) 
print(dim(edge_df_final)) # dim = (28,7)

# COMMAND ----------

# Creating a weighted graph
g_filtered = graph_from_data_frame(data.frame(
  from = edge_df_final$start,
  to = edge_df_final$end,
  weight = edge_df_final$n_patient), directed = TRUE)
print(length(V(g_filtered)$name)) # 17
# Get all paths from all vertices to vertex G30 (remove edges that cannot lead to G30 b/c we only care about paths to G30)
allPathsTo_G30 = lapply(1:length(V(g_filtered)), function(v) {
  if(V(g_filtered)$name[v] != "G30") { # Ensure not to find paths from G30 to G30
    paths = all_simple_paths(g_filtered, from = V(g_filtered)$name[v], to = "G30")
    return(lapply(paths, function(p) V(g_filtered)[p]$name))
  }
})
allPathsTo_G30 = Filter(function(x) !is.null(x) && length(x) > 0, allPathsTo_G30)
allPathsTo_G30 = unlist(allPathsTo_G30, recursive = FALSE)
# Find the maximum length of the paths
max_length = max(sapply(allPathsTo_G30, length))
# Pad the shorter paths with NAs to make all paths equal in length
allPathsTo_G30 = lapply(allPathsTo_G30, function(x) {
  length(x) = max_length  
  return(x)
})
# Convert the list of paths to a data frame
paths_df = as.data.frame(do.call(rbind, allPathsTo_G30), stringsAsFactors = FALSE) %>% filter(!is.na(V3))
print(dim(paths_df)) # dim = (71,7)

# COMMAND ----------

# restrict to nodes to G30
node_to_G30 = unique(c(paths_df$V1, paths_df$V2, paths_df$V3, paths_df$V4, paths_df$V5, paths_df$V6, paths_df$V7))
node_to_G30 = node_to_G30[!is.na(node_to_G30)]
# summarize results
filter_edges = as.data.frame(get.edgelist(g)) %>% rename("start" = "V1", "end" = "V2") %>%
  inner_join(pat_start_end_summary) %>% 
  filter(start %in% node_to_G30 & end %in% node_to_G30)

# COMMAND ----------

filter_edges %>% arrange(n_patient) %>% head(20)

# COMMAND ----------

filter_edges %>% filter(end == "G93") %>% arrange(desc(n_patient)) %>% head()

# COMMAND ----------

download_path = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/output/2024_09_16/"
# save(filter_edges, file = paste0(download_path, "filter_edges_", cluster_i, ".rda"))
write.csv(filter_edges, file = paste0(download_path, "filter_edges_", cluster_i, ".csv"), row.names = FALSE)

# COMMAND ----------

save(filter_edges, file = paste0(raw_data_path, "clusters/filter_edges_", cluster_i, ".rda"))

# COMMAND ----------

# Plot network backbone
plot_df = filter_edges
nodes_traj = as_tibble(union(names(table(plot_df$start)), names(table(plot_df$end)))) %>% 
  rowid_to_column("id") %>% rename(label = value) %>% 
  mutate(class = substr(label, start = 1, stop = 1))
per_route_traj = plot_df %>% select(start, end, n_patient) %>% rename(weight = n_patient)
edges_traj = per_route_traj %>% 
  left_join(nodes_traj, by = c("start" = "label")) %>% rename(from = id) %>% 
  left_join(nodes_traj, by = c("end" = "label")) %>% rename(to = id) %>% select(from, to, weight)
routes_tidy_traj = tbl_graph(nodes = nodes_traj, edges = edges_traj, directed = TRUE) %>% 
  activate(edges) %>% arrange(desc(weight))
# plot
set.seed(1234)
# pdf(paste0("/Workspace/Users/mingzhoufu@mednet.ucla.edu/outputs/cluster/", cluster_i, "_backbone.pdf") , 
#     width = 6, height = 4)
# ggraph(routes_tidy_traj, layout = 'graphopt') + 
#   geom_node_point(aes(colour = class), size = 7) +
#   geom_edge_link(aes(width = weight), alpha = 0.8, arrow = arrow(length = unit(2, "mm")), end_cap = circle(1, "mm")) + 
#   scale_edge_width(range = c(0.2, 2)) + geom_node_text(aes(label = label), repel = TRUE) +
#   labs(edge_width = "Number of cases") + theme_graph(base_family = 'Helvetica')
# dev.off()
ggraph(routes_tidy_traj, layout = 'graphopt') + 
  geom_node_point(aes(colour = class), size = 14) +
  geom_edge_link(aes(width = weight), alpha = 0.8, arrow = arrow(length = unit(2, "mm")), end_cap = circle(2, "mm")) + 
  scale_edge_width(range = c(0.2, 2)) + geom_node_text(aes(label = label), repel = TRUE, size = 6) +
  labs(edge_width = "Number of cases") + theme_graph(base_family = 'Calibri')

# COMMAND ----------

# MAGIC %md
# MAGIC #### Central nodes

# COMMAND ----------

central_nodes = combined_modularity_df %>% filter(desc %in% final_node & desc != "G30") %>% arrange(desc(vitality)) %>% head(10) %>% pull(desc)
# Add descriptions
icd10cm2016_3digit = icd10cm2016 %>% filter(nchar(code) == 3) %>% 
  mutate(code = as.character(code))
central_nodes_df = cbind(c(1:5), central_nodes) %>% as.data.frame() %>% 
  left_join(icd10cm2016_3digit, by = c("central_nodes" = "three_digit")) %>% 
  select(central_nodes, long_desc, sub_chapter, chapter)
central_nodes_df

# COMMAND ----------

# MAGIC %md
# MAGIC #### Trajectories in backbones

# COMMAND ----------

# We will use their full trajectories (without greedy concatenation, i.e., subsets of their full trajectories)
id_lst = unique(pat_cluster_i_long$PatientID)
print(length(id_lst)) # 3223
# combine full results
folder_path = paste0(raw_data_path, "trajectory_clean/no_impute/short_traj/pat_full/") 
# List all .rda files from the folder
files = list.files(path = folder_path, pattern = "\\.rda$")
print(length(files)) # 4191

# COMMAND ----------

num_files = sapply(strsplit(files, "_"), `[`, 1) %>% as.numeric() %>% sort()
to_load_files = c()
for (id in id_lst) {
  # Use grep to find elements containing patient_id
  matching_file = files[grep(paste0("_", id, ".rda"), files)]
  to_load_files = c(to_load_files, matching_file) %>% unique()
}
print(length(to_load_files)) # should be exact same as length(id_lst)
# load in their full files
for (i in 1:length(to_load_files)){
  print(i)
  load(paste0(folder_path, to_load_files[i]))
  paths_pat_full_final$PatientID = as.character(paths_pat_full_final$PatientID)
  print(dim(paths_pat_full_final))
  if (i == 1) {
    paths_pat_combine_full = paths_pat_full_final
  } else {
    paths_pat_combine_full = bind_rows(paths_pat_combine_full, paths_pat_full_final) %>% unique()
  }
}
print(dim(paths_pat_combine_full)) # dim = (450308,17)
names(impute_wide_3step)= c("PatientID", "V1", "V2", "V3")
paths_pat_combine_cluster = bind_rows(paths_pat_combine_full, impute_wide_3step)
print(dim(paths_pat_combine_cluster)) # dim = (456899,17)
save(paths_pat_combine_cluster, file = paste0(raw_data_path, "trajectory_clean/final/paths_pat_combine_cluster_", cluster_i, ".rda"))

# COMMAND ----------

load(file = paste0(raw_data_path, "trajectory_clean/final/paths_pat_combine_cluster_", cluster_i, ".rda"))
pat_cluster_i_wide_full = paths_pat_combine_cluster %>% filter(!is.na(V3)) %>% unique() %>% 
  group_by(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15) %>% 
  summarise(n_patient = length(unique(PatientID))) %>% ungroup() 
pat_cluster_i_wide_full_short = pat_cluster_i_wide_full %>% filter(is.na(V8)) %>%
  select(V1, V2, V3, V4, V5, V6, V7, n_patient)
final_pat_count = pat_cluster_i_wide_full_short %>% inner_join(paths_df) %>% arrange(desc(n_patient))
head(final_pat_count, 20)

# COMMAND ----------

final_pat_count %>% filter(V2 == "G93" | V3 == "G93")

# COMMAND ----------

# MAGIC %md
# MAGIC #### Look at >3-digit ICDs

# COMMAND ----------

length(id_lst)

# COMMAND ----------

test_icd = "G93"
pat_detail_diag = AD_cluster_ICD_filter %>% filter(PatientID %in% id_lst & ICD_3digit == test_icd) %>% 
  group_by(ICD) %>% summarise(freq = length(unique(PatientID))) 
pie_desc = pat_detail_diag %>% 
  # remove "." in the code
  mutate(code = str_remove(ICD, "\\.")) %>% left_join(icd10cm2016_short, by = c("code" = "code")) %>% 
  drop_na() %>% select(-code) %>% filter(ICD != test_icd) %>% mutate(prop = round(freq / length(id_lst)*100, 2)) %>%
  mutate(freq_prop = paste0(freq, " (", prop, "%)")) %>% select(ICD, long_desc, freq_prop)
pie_desc

# COMMAND ----------

tt = AD_cluster_ICD_filter %>% filter(PatientID %in% id_lst & ICD_3digit == test_icd) %>% pull(PatientID) %>% unique() %>% length()
pat_detail_diag$freq/tt*100

# COMMAND ----------

tt

# COMMAND ----------

pat_detail_diag

# COMMAND ----------

# MAGIC %md
# MAGIC ## Cluster 4: Neurodegeneration trajectory

# COMMAND ----------

cluster_i = 4

# COMMAND ----------

# MAGIC %md
# MAGIC #### Raw network

# COMMAND ----------

pat_cluster_i_long = pat_traj_final_cluster %>% filter(Kmeans == cluster_i) %>% 
  select(PatientID, paste0("step", c(1:9))) %>% 
  pivot_longer(cols = paste0("step", c(1:9)), names_to = "step", values_to = "ICD_3digit") %>% drop_na() %>% 
  mutate(step = as.numeric(str_sub(step, start = 5, end = -1L))) 
id_lst = unique(pat_cluster_i_long$PatientID)
print(length(id_lst)) # 1502
for (n in 1:length(id_lst)) {
  patient_id = id_lst[n]
  # print(paste0("Patinet ", n, ": ", patient_id))
  pat_rec = pat_cluster_i_long %>% filter(PatientID == patient_id) %>% 
    group_by(step) %>% mutate(traj = row_number()) %>% ungroup()
  n_traj = max(pat_rec$traj)
  for (i in 1:n_traj) {
    icd_vec = pat_rec %>% filter(traj == i) %>% pull(ICD_3digit)
    start_end_sub = data.frame(start = icd_vec[-length(icd_vec)], end = icd_vec[-1]) %>% 
      mutate(PatientID = patient_id, traj = i, n = 1)
    if (i == 1) {patient_start_end = start_end_sub}
    else {patient_start_end = rbind(patient_start_end, start_end_sub)}
  }
  if (n == 1) {start_end_full = patient_start_end}
  else {start_end_full = rbind(start_end_full, patient_start_end)}
}
pat_start_end_summary = start_end_full %>% select(-traj) %>% unique() %>% 
  group_by(start, end) %>% summarise(n_patient = n()) %>% filter(start != "G30") %>% ungroup()
print(dim(pat_start_end_summary)) # dim = (1993,3)

# COMMAND ----------

pat_cluster_i_long %>% filter(ICD_3digit == "G31") %>% pull(PatientID) %>% unique() %>% length()

# COMMAND ----------

# filter based on significant pairs and #pat followed (>=0.5%)
pat_start_end_summary_sig = pat_start_end_summary %>% 
  left_join(risk_pair_sig %>% select(exposure_ICD, outcome_ICD, coef) %>% 
  rename(start = exposure_ICD, end = outcome_ICD), by = c("start", "end")) %>% 
  filter(!is.na(coef)) %>% select(-coef) %>% filter(n_patient >= length(id_lst)*0.005)
print(dim(pat_start_end_summary_sig)) # dim = (112,3)
# Creating a weighted graph
g = graph_from_data_frame(data.frame(
  from = pat_start_end_summary_sig$start,
  to = pat_start_end_summary_sig$end,
  weight = pat_start_end_summary_sig$n_patient), directed = TRUE)
print(length(V(g)$name)) # 51
# Get all paths from all vertices to vertex G30 (remove edges that cannot lead to G30 b/c we only care about paths to G30)
allPathsTo_G30 = lapply(1:length(V(g)), function(v) {
  if(V(g)$name[v] != "G30") { # Ensure not to find paths from G30 to G30
    paths = all_simple_paths(g, from = V(g)$name[v], to = "G30")
    return(lapply(paths, function(p) V(g)[p]$name))
  }
})
allPathsTo_G30 = Filter(function(x) !is.null(x) && length(x) > 0, allPathsTo_G30)
allPathsTo_G30 = unlist(allPathsTo_G30, recursive = FALSE)
# Find the maximum length of the paths
max_length = max(sapply(allPathsTo_G30, length))
# Pad the shorter paths with NAs to make all paths equal in length
allPathsTo_G30 = lapply(allPathsTo_G30, function(x) {
  length(x) = max_length  
  return(x)
})
# Convert the list of paths to a data frame
paths_df = as.data.frame(do.call(rbind, allPathsTo_G30), stringsAsFactors = FALSE) %>% filter(!is.na(V3))
print(dim(paths_df)) # dim = (4242,12)

# COMMAND ----------

# restrict to nodes to G30
node_to_G30 = unique(c(paths_df$V1, paths_df$V2, paths_df$V3, paths_df$V4, paths_df$V5, paths_df$V6, 
  paths_df$V7, paths_df$V8, paths_df$V9, paths_df$V10, paths_df$V11, paths_df$V12))
node_to_G30 = node_to_G30[!is.na(node_to_G30)]
# summarize results
filter_edges = as.data.frame(get.edgelist(g)) %>% rename("start" = "V1", "end" = "V2") %>%
  inner_join(pat_start_end_summary) %>% 
  filter(start %in% node_to_G30 & end %in% node_to_G30)
# Plot network backbone
plot_df = filter_edges
nodes_traj = as_tibble(union(names(table(plot_df$start)), names(table(plot_df$end)))) %>% 
  rowid_to_column("id") %>% rename(label = value) %>% 
  mutate(class = substr(label, start = 1, stop = 1))
per_route_traj = plot_df %>% select(start, end, n_patient) %>% rename(weight = n_patient)
edges_traj = per_route_traj %>% 
  left_join(nodes_traj, by = c("start" = "label")) %>% rename(from = id) %>% 
  left_join(nodes_traj, by = c("end" = "label")) %>% rename(to = id) %>% select(from, to, weight)
routes_tidy_traj = tbl_graph(nodes = nodes_traj, edges = edges_traj, directed = TRUE) %>% 
  activate(edges) %>% arrange(desc(weight))

# COMMAND ----------

# plot
set.seed(1234)
# pdf(paste0("/Workspace/Users/mingzhoufu@mednet.ucla.edu/outputs/cluster/", cluster_i, "_backbone.pdf") , 
#     width = 6, height = 4)
ggraph(routes_tidy_traj, layout = 'graphopt') + 
  geom_node_point(aes(colour = class), size = 7) +
  geom_edge_link(aes(width = weight), alpha = 0.8, arrow = arrow(length = unit(2, "mm")), end_cap = circle(1, "mm")) + 
  scale_edge_width(range = c(0.2, 2)) + geom_node_text(aes(label = label), repel = TRUE) +
  labs(edge_width = "Number of cases") + theme_graph(base_family = 'Helvetica')
# dev.off()

# COMMAND ----------

# MAGIC %md
# MAGIC ### Backbone extraction - Modularity Vitality

# COMMAND ----------

# start from significant pairs 
pat_start_end_summary_sig = pat_start_end_summary %>% 
  left_join(risk_pair_sig %>% select(exposure_ICD, outcome_ICD, coef) %>% 
  rename(start = exposure_ICD, end = outcome_ICD), by = c("start", "end")) %>% 
  filter(!is.na(coef)) %>% select(-coef) %>% filter(n_patient >= length(id_lst)*0.005)
print(dim(pat_start_end_summary_sig)) # dim = (112,3)
# Creating a weighted graph
g = graph_from_data_frame(data.frame(
  from = pat_start_end_summary_sig$start,
  to = pat_start_end_summary_sig$end,
  weight = pat_start_end_summary_sig$n_patient), directed = TRUE)
print(length(V(g)$name)) # 51

# COMMAND ----------

# MAGIC %md
# MAGIC #### Sort modularity

# COMMAND ----------

# Initialize lists to store modularity vitality scores for nodes
node_vitality = numeric(vcount(g))
# Compute modularity vitality for each node
for (v in V(g)) {
  node_vitality[v] = modularity_vitality_node(g, v)
}
# create a vitality dataframe (combine nodes and edges)
node_df = data.frame(type = "node", id = V(g), vitality = abs(node_vitality)) %>% as.data.frame() %>%
  mutate(vitality = as.numeric(vitality), id = as.numeric(id))
rownames(node_df) = NULL
node_desc = data.frame(id = V(g), desc = V(g)$name) %>% as.data.frame() %>%
  mutate(id = as.numeric(id), desc = as.character(desc))
node_df = node_df %>% left_join(node_desc) %>% arrange(vitality)
print(dim(node_df)) # dim = (51,4)
tail(node_df)

# COMMAND ----------

# Initialize lists to store modularity vitality scores for edges
edge_vitality = numeric(ecount(g))
# Compute modularity vitality for each edge
for (e in E(g)) {
  edge_vitality[e] = modularity_vitality_edge(g, e)
}
edge_df = data.frame(type = "edge", id = E(g), vitality = abs(edge_vitality)) %>% as.data.frame() %>%
  mutate(vitality = as.numeric(vitality), id = as.numeric(id))
edge_list = get.edgelist(g, names = TRUE)
edge_pairs = apply(edge_list, 1, function(x) paste(x[1], x[2], sep = " -> "))
edge_desc = data.frame(desc = edge_pairs, stringsAsFactors = FALSE) %>% rownames_to_column(var = "id") %>%
  mutate(id = as.numeric(id), desc = as.character(desc))
edge_df = edge_df %>% left_join(edge_desc) %>% arrange(vitality) %>%
  separate(desc, into = c("start", "end"), sep = " -> ", remove = FALSE)
print(dim(edge_df)) # dim = (112,6)
print(tail(edge_df))

# COMMAND ----------

edge_df_short = edge_df %>% select(-c(start, end))
combined_modularity_df = rbind(node_df, edge_df_short) %>% as.data.frame() %>% arrange(vitality)
print(dim(combined_modularity_df)) # dim = (163,4)
print(tail(combined_modularity_df))

# COMMAND ----------

# check full modularity change with one node/edge removed at a time
g_rep = g
modularity_lst = compute_modularity(g_rep)
for (i in 1:162) {
  element_to_remove = combined_modularity_df[i,]
  if (element_to_remove$type == "node") {
    node_to_remove = element_to_remove$desc
    g_rep = delete_vertices(g_rep, node_to_remove)
  } else if (element_to_remove$type == "edge") {
    edge_to_remove = element_to_remove$desc
    edge_list = get.edgelist(g_rep, names = TRUE)
    edge_pairs = apply(edge_list, 1, function(x) paste(x[1], x[2], sep = " -> "))
    edge_desc = data.frame(desc = edge_pairs, stringsAsFactors = FALSE) %>% rownames_to_column(var = "id") %>%
      mutate(id = as.numeric(id), desc = as.character(desc))
    index_edge_rmv = which(edge_desc$desc == edge_to_remove)
    g_rep = delete_edges(g_rep, index_edge_rmv)
  }
  modularity_lst = c(modularity_lst, compute_modularity(g_rep))
}
mod_change_df = cbind(c(0:162), modularity_lst) %>% as.data.frame()
names(mod_change_df) = c("threshold", "modularity")
plot(mod_change_df, type = "o", col = "blue", xlab = "Step", ylab = "Modularity")

# COMMAND ----------

mod_change_df %>% filter(threshold < 120) %>% arrange(desc(modularity)) %>% head() # N = 101

# COMMAND ----------

node_edge_filter = combined_modularity_df[101:163,]
print(dim(node_edge_filter)) # dim = (63,4)
print(tail(node_edge_filter))
node_edge_select = node_edge_filter %>% pull(desc)
print(table(node_edge_filter$type))
edge_df_filter = edge_df %>% filter(desc %in% node_edge_filter$desc)
edge_filter_node = unique(c(edge_df_filter$start, edge_df_filter$end))
node_filter = node_edge_filter %>% filter(type == "node") %>% pull(desc) %>% unique()
final_node = c("G30", intersect(node_filter, edge_filter_node)) %>% unique()
print(length(final_node)) # 24
edge_df_final = edge_df_filter %>% filter(start %in% final_node & end %in% final_node) %>% left_join(pat_start_end_summary_sig) 
print(dim(edge_df_final)) # dim = (25,7)

# COMMAND ----------

# Creating a weighted graph
g_filtered = graph_from_data_frame(data.frame(
  from = edge_df_final$start,
  to = edge_df_final$end,
  weight = edge_df_final$n_patient), directed = TRUE)
print(length(V(g_filtered)$name)) # 22
# Get all paths from all vertices to vertex G30 (remove edges that cannot lead to G30 b/c we only care about paths to G30)
allPathsTo_G30 = lapply(1:length(V(g_filtered)), function(v) {
  if(V(g_filtered)$name[v] != "G30") { # Ensure not to find paths from G30 to G30
    paths = all_simple_paths(g_filtered, from = V(g_filtered)$name[v], to = "G30")
    return(lapply(paths, function(p) V(g_filtered)[p]$name))
  }
})
allPathsTo_G30 = Filter(function(x) !is.null(x) && length(x) > 0, allPathsTo_G30)
allPathsTo_G30 = unlist(allPathsTo_G30, recursive = FALSE)
# Find the maximum length of the paths
max_length = max(sapply(allPathsTo_G30, length))
# Pad the shorter paths with NAs to make all paths equal in length
allPathsTo_G30 = lapply(allPathsTo_G30, function(x) {
  length(x) = max_length  
  return(x)
})
# Convert the list of paths to a data frame
paths_df = as.data.frame(do.call(rbind, allPathsTo_G30), stringsAsFactors = FALSE) %>% filter(!is.na(V3))
print(dim(paths_df)) # dim = (26,5)

# COMMAND ----------

# restrict to nodes to G30
node_to_G30 = unique(c(paths_df$V1, paths_df$V2, paths_df$V3, paths_df$V4, paths_df$V5))
node_to_G30 = node_to_G30[!is.na(node_to_G30)]
# summarize results
filter_edges = as.data.frame(get.edgelist(g)) %>% rename("start" = "V1", "end" = "V2") %>%
  inner_join(pat_start_end_summary) %>% 
  filter(start %in% node_to_G30 & end %in% node_to_G30)

# COMMAND ----------

filter_edges %>% arrange(n_patient) %>% head(10)

# COMMAND ----------

filter_edges %>% filter(end == "G31") %>% arrange(desc(n_patient)) %>% head()

# COMMAND ----------

download_path = "/Workspace/Users/mingzhoufu@mednet.ucla.edu/output/2024_09_16/"
# save(filter_edges, file = paste0(download_path, "filter_edges_", cluster_i, ".rda"))
write.csv(filter_edges, file = paste0(download_path, "filter_edges_", cluster_i, ".csv"), row.names = FALSE)

# COMMAND ----------

save(filter_edges, file = paste0(raw_data_path, "clusters/filter_edges_", cluster_i, ".rda"))

# COMMAND ----------

# Plot network backbone
plot_df = filter_edges
nodes_traj = as_tibble(union(names(table(plot_df$start)), names(table(plot_df$end)))) %>% 
  rowid_to_column("id") %>% rename(label = value) %>% 
  mutate(class = substr(label, start = 1, stop = 1))
per_route_traj = plot_df %>% select(start, end, n_patient) %>% rename(weight = n_patient)
edges_traj = per_route_traj %>% 
  left_join(nodes_traj, by = c("start" = "label")) %>% rename(from = id) %>% 
  left_join(nodes_traj, by = c("end" = "label")) %>% rename(to = id) %>% select(from, to, weight)
routes_tidy_traj = tbl_graph(nodes = nodes_traj, edges = edges_traj, directed = TRUE) %>% 
  activate(edges) %>% arrange(desc(weight))
# plot
set.seed(1234)
# pdf(paste0("/Workspace/Users/mingzhoufu@mednet.ucla.edu/outputs/cluster/", cluster_i, "_backbone.pdf") , 
#     width = 6, height = 4)
# ggraph(routes_tidy_traj, layout = 'graphopt') + 
#   geom_node_point(aes(colour = class), size = 7) +
#   geom_edge_link(aes(width = weight), alpha = 0.8, arrow = arrow(length = unit(2, "mm")), end_cap = circle(1, "mm")) + 
#   scale_edge_width(range = c(0.2, 2)) + geom_node_text(aes(label = label), repel = TRUE) +
#   labs(edge_width = "Number of cases") + theme_graph(base_family = 'Helvetica')
# dev.off()
ggraph(routes_tidy_traj, layout = 'graphopt') + 
  geom_node_point(aes(colour = class), size = 14) +
  geom_edge_link(aes(width = weight), alpha = 0.8, arrow = arrow(length = unit(2, "mm")), end_cap = circle(2, "mm")) + 
  scale_edge_width(range = c(0.2, 2)) + geom_node_text(aes(label = label), repel = TRUE, size = 6) +
  labs(edge_width = "Number of cases") + theme_graph(base_family = 'Calibri')

# COMMAND ----------

# MAGIC %md
# MAGIC #### Central nodes

# COMMAND ----------

central_nodes = combined_modularity_df %>% filter(desc %in% final_node & desc != "G30") %>% arrange(desc(vitality)) %>% head(10) %>% pull(desc)
# Add descriptions
icd10cm2016_3digit = icd10cm2016 %>% filter(nchar(code) == 3) %>% 
  mutate(code = as.character(code))
central_nodes_df = cbind(c(1:5), central_nodes) %>% as.data.frame() %>% 
  left_join(icd10cm2016_3digit, by = c("central_nodes" = "three_digit")) %>% 
  select(central_nodes, long_desc, sub_chapter, chapter)
central_nodes_df

# COMMAND ----------

# MAGIC %md
# MAGIC #### Trajectories in backbones

# COMMAND ----------

# We will use their full trajectories (without greedy concatenation, i.e., subsets of their full trajectories)
id_lst = unique(pat_cluster_i_long$PatientID)
print(length(id_lst)) # 1502
# combine full results
folder_path = paste0(raw_data_path, "trajectory_clean/no_impute/short_traj/pat_full/") 
# List all .rda files from the folder
files = list.files(path = folder_path, pattern = "\\.rda$")
print(length(files)) # 4191

# COMMAND ----------

num_files = sapply(strsplit(files, "_"), `[`, 1) %>% as.numeric() %>% sort()
to_load_files = c()
for (id in id_lst) {
  # Use grep to find elements containing patient_id
  matching_file = files[grep(paste0("_", id, ".rda"), files)]
  to_load_files = c(to_load_files, matching_file) %>% unique()
}
print(length(to_load_files)) # should be exact same as length(id_lst)
# load in their full files
for (i in 1:length(to_load_files)){
  print(i)
  load(paste0(folder_path, to_load_files[i]))
  paths_pat_full_final$PatientID = as.character(paths_pat_full_final$PatientID)
  print(dim(paths_pat_full_final))
  if (i == 1) {
    paths_pat_combine_full = paths_pat_full_final
  } else {
    paths_pat_combine_full = bind_rows(paths_pat_combine_full, paths_pat_full_final) %>% unique()
  }
}
print(dim(paths_pat_combine_full)) # dim = (467942,16)
names(impute_wide_3step)= c("PatientID", "V1", "V2", "V3")
paths_pat_combine_cluster = bind_rows(paths_pat_combine_full, impute_wide_3step)
print(dim(paths_pat_combine_cluster)) # dim = (474533,16)
save(paths_pat_combine_cluster, file = paste0(raw_data_path, "trajectory_clean/final/paths_pat_combine_cluster_", cluster_i, ".rda"))

# COMMAND ----------

# load(file = paste0(raw_data_path, "trajectory_clean/final/paths_pat_combine_cluster_", cluster_i, ".rda"))
pat_cluster_i_wide_full = paths_pat_combine_cluster %>% filter(!is.na(V3)) %>% unique() %>% 
  group_by(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15) %>% 
  summarise(n_patient = length(unique(PatientID))) %>% ungroup() 
pat_cluster_i_wide_full_short = pat_cluster_i_wide_full %>% filter(is.na(V6)) %>%
  select(V1, V2, V3, V4, V5, n_patient)
final_pat_count = pat_cluster_i_wide_full_short %>% inner_join(paths_df) %>% arrange(desc(n_patient))
head(final_pat_count, 20)

# COMMAND ----------

# MAGIC %md
# MAGIC #### Look at >3-digit ICDs

# COMMAND ----------

length(id_lst)

# COMMAND ----------

test_icd = "G31"
pat_detail_diag = AD_cluster_ICD_filter %>% filter(PatientID %in% id_lst & ICD_3digit == test_icd) %>% 
  group_by(ICD) %>% summarise(freq = length(unique(PatientID))) 
pie_desc = pat_detail_diag %>% 
  # remove "." in the code
  mutate(code = str_remove(ICD, "\\.")) %>% left_join(icd10cm2016_short, by = c("code" = "code")) %>% 
  drop_na() %>% select(-code) %>% filter(ICD != test_icd) %>% mutate(prop = round(freq / length(id_lst)*100, 2)) %>%
  mutate(freq_prop = paste0(freq, " (", prop, "%)")) %>% select(ICD, long_desc, freq_prop)
pie_desc

# COMMAND ----------

AD_cluster_ICD_filter %>% filter(PatientID %in% id_lst & ICD_3digit == test_icd) %>% pull(PatientID) %>% unique() %>% length()

# COMMAND ----------

pat_detail_diag

# COMMAND ----------

# MAGIC %md
# MAGIC # New Added

# COMMAND ----------

load(file = paste0(raw_data_path, "FG_test/other_risk/diag_beforeAD_filter_0610.rda"))
print(dim(diag_beforeAD_filter)) # dim = (393601,4)
head(diag_beforeAD_filter)

# COMMAND ----------

cluster_i = 4

# COMMAND ----------

pat_cluster_i_long = pat_traj_final_cluster %>% filter(Kmeans == cluster_i) %>% 
  select(PatientID, paste0("step", c(1:9))) %>% 
  pivot_longer(cols = paste0("step", c(1:9)), names_to = "step", values_to = "ICD_3digit") %>% drop_na() %>% 
  mutate(step = as.numeric(str_sub(step, start = 5, end = -1L))) 
id_lst = unique(pat_cluster_i_long$PatientID)
print(length(id_lst)) # 1448

# COMMAND ----------

icd_1 = "G31"
icd_2 = "F02"

pat_record = diag_beforeAD_filter %>% filter(PatientID %in% id_lst) %>%
  filter(ICD_3digit %in% c(icd_1, icd_2)) %>%
  pivot_wider(names_from = ICD_3digit, values_from = EncounterDate, values_fill = NA_Date_) %>% 
  select(-AD_Date) %>% mutate(time = as.numeric(get(icd_2) - get(icd_1))/365.25) %>%
  drop_na() %>% filter(time >= 0)
print(dim(pat_record))
print(head(pat_record))
print(mean(pat_record$time, na.rm=T))

# COMMAND ----------

