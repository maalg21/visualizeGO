## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

## ----install, eval=FALSE------------------------------------------------------
# # install.packages("devtools")
# devtools::install_github("maalg21/visualizeGO", force = TRUE)

## ----load_data----------------------------------------------------------------
library(visualizeGO)
library(dplyr)
library(igraph)

data <- as.data.frame(readxl::read_xlsx(system.file("extdata", "GOTerms.xlsx", package = "visualizeGO")))
data <- data[,c("Category", "ID", "Padj")]
GO_BP <- data[data$Category == "BP",]

## ----calculate_similarity-----------------------------------------------------
similarity_matrix <- calculate_wang(input_terms = GO_BP$ID,  ontology = "BP", # We select the category to which our GO-terms belong.
OrgDb = "org.Hs.eg.db") # We use human annotation as a reference.

## ----filter_terms-------------------------------------------------------------
filtered_terms <- filter_terms(similarity_matrix = similarity_matrix, threshold = 0.1)
similarity_matrix2 <- filtered_terms$similarity_matrix

## ----determine_nbclusters-----------------------------------------------------
set.seed(1234) # For reproducibility
# For the similarity matrix NOT filtered:
determine_nbclusters(similarity_matrix = similarity_matrix)
# Based on the Elbow method, the optimal number of clusters (k) is: 20 
# Based on the Silhouette method, the optimal number of clusters (k) is: 54

## ----determine_nbclusters2----------------------------------------------------
# For the filtered similarity matrix:
determine_nbclusters(similarity_matrix = similarity_matrix2)
# Based on the Elbow method, the optimal number of clusters (k) is: 2 
# Based on the Silhouette method, the optimal number of clusters (k) is: 34

## ----cluster------------------------------------------------------------------
cluster <- clusterGO(similarity_matrix = similarity_matrix, nb_clusters = 54)

## ----cluster2-----------------------------------------------------------------
cluster2 <- clusterGO(similarity_matrix = similarity_matrix2, nb_clusters = 34)

## ----table--------------------------------------------------------------------
colors <- generate_pastel_colors(n = 54) # This function was only created to generate a list of pastel colours of the number we determine 😊
generate_cluster_table(cluster_output = cluster, col_palette = colors,  text_color = "black", file_name = NULL) # If you specify a name for the file, it will be saved as a PNG.

## ----table2-------------------------------------------------------------------
colors <- generate_pastel_colors(n = 34)
generate_cluster_table(cluster_output = cluster2, col_palette = colors, 
text_color = "black", file_name = NULL) # If you specify a name for the file, it will be saved as a PNG.

## ----scater_plot--------------------------------------------------------------
scatterGO(similarity_matrix = similarity_matrix, cluster = cluster, colors = generate_pastel_colors(n = 54), title = "Distance Between GO-Terms", labels = T, size = "padj", scores = setNames(-log10(GO_BP$Padj), GO_BP$ID))

## ----treemap_plot-------------------------------------------------------------
treeMap(cluster, size = "padj", scores = setNames(-log10(GO_BP$Padj), GO_BP$ID), title = "Distance Between GO-Terms", colors = generate_pastel_colors(n = 54))

## ----compare_GOlists----------------------------------------------------------
# This is the second list of GO-terms
data2 <- as.data.frame(readxl::read_xlsx(system.file("extdata", "GOTerms2.xlsx", package = "visualizeGO")))

# We filter 10 GO-terms for each list, to make it more easy to understand
GO_BP1 <- data[data$Category == "BP", c("Category", "ID")] %>% top_n(n = 10)
GO_BP2 <- data2[data2$Category == "BP",] %>% top_n(n = 10)

compareGO(comparison = "GO", list1 = GO_BP1$ID, list2 = GO_BP2$ID, ontology = "BP", OrgDb = "org.Hs.eg.db", method = "Wang", plot = T, low = "white", high = "red3", labs = c("GO List 1", "GO List 2"), cex = 3, cex_axis = 10)

## ----compare_clusters---------------------------------------------------------
compareGO(comparison = "cluster", list1 = cluster$clusters, list2 = cluster2$clusters, ontology = "BP", OrgDb = "org.Hs.eg.db", method = "Wang", combine = "BMA", plot = T, low = "white", high = "red3", labs = NULL, cex = 3, cex_axis = 10)

## ----hierarchical_graph-------------------------------------------------------
graph <- familyGO(cluster, go_sim_object = NULL)

## ----visualizeGO_plot---------------------------------------------------------
visualizeGO(cluster = cluster, graph = graph, shape = "circle", min_node_size = 1, max_node_size = 10, layout = "tree", title = "Alonso-García et al. (2023)", labs = "Input Terms", col_palette = generate_pastel_colors(54), verbose = "some", legend = T, ID = F)

## ----visualizeGO_plot_filtered------------------------------------------------
visualizeGO(cluster = cluster, graph = graph,
selected_cluster = c(13,16,18,19), 
shape = "circle", min_node_size = 2.5, max_node_size = 10, 
layout = "tree", col_palette = colors[c(13,16,18,19)], 
title = "Adipose Tissue Development related clusters",
verbose = "some", legend = T, ID = T,
labs = "Input Terms")

