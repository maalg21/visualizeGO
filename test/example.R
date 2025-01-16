devtools::install_github("maalg21/visualizeGO", force = T)
library(visualizeGO)

# Used data ----
data <- as.data.frame(readxl::read_xlsx(system.file("extdata", "GOTerms.xlsx",
                                                    package = "visualizeGO")))
data <- data[,c("Category", "ID")]
GO_BP <- data[data$Category == "BP",]

# Building the original graph ----
graph <- build_hierarchical_graph(go_list1 = GO_BP$ID, nb_lists = "single",
                                  go_sim_object = NULL)
expanded_graph <- expand_graph(graph = graph,
                               go_list1 = GO_BP$ID,
                               nb_lists = "single")
final_graph <- retain_ancestors_above_input_terms(graph = expanded_graph,
                                                  go_list1 = GO_BP$ID,
                                                  nb_lists = "single")

# Semantic similarity clustering ----
similarity_matrix <- calculate_wang(graph = final_graph,
                                    ontology = "BP",
                                    orgdb = "org.Hs.eg.db")
cluster <- cluster_go_terms(method_type = "similarity", method = "wang",
                            orgdb = "org.Hs.eg.db", ontology = "BP", similarity_matrix = similarity_matrix,
                            graph = final_graph, nb_clusters = NULL, k_range = 2:10)

# Whole plot ----
colors <- generate_pastel_colors(n = 7)
visualize_go_hierarchy(go_list1 = GO_BP$ID,
                       go_list2 = NULL,
                       nb_lists = "single",
                       go_sim_object = NULL,
                       shape1 = "circle", # By default, it would be circles
                       shape2 = NULL,
                       ontology = "BP",
                       simplification = T,
                       min_node_size = 1, max_node_size = 10, # This is optional and arbitrary
                       layout = "tree",
                       clustering = T, clusters = cluster,
                       col_palette = colors,
                       legend = T, labs = "GO-Terms",
                       verbose = "some", # Just to know more about the process that is ocurring
                       save_plot = T, # For saving the following plot
                       PNG = "plot1.png")

# Filtering to obtain only nodes of Cluster 5 ----
filter_and_visualize_cluster(clusters = cluster,
                             selected_cluster = 5,
                             ontology = "BP",
                             layout = "tree", col_palette = "#D9C9FF",
                             min_node_size = 1, max_node_size = 10,
                             save_plot = T, PNG = "plot2.png",
                             verbose = "some", legend = T)

# # Filtering to obtain nodes of Cluster 6 & 7 ----
filter_and_visualize_cluster(clusters = cluster,
                             selected_cluster = c(6,7),
                             ontology = "BP",
                             layout = "tree",
                             col_palette = generate_pastel_colors(n = 7)[c(6,7)],
                             min_node_size = 1, max_node_size = 10,
                             save_plot = T, PNG = "plot3.png",
                             verbose = "some", legend = T)

# Comparing two lists of GO-terms ----
data2 <- as.data.frame(readxl::read_xlsx(system.file("extdata", "GOTerms2.xlsx",
                                                     package = "visualizeGO")))

GO_BP1 <- data[data$Category == "BP",] %>% top_n(n = 10)
GO_BP2 <- data2[data2$Category == "BP",] %>% top_n(n = 10)

go_similarity_heatmap(go_list1 = GO_BP1$ID, go_list2 = GO_BP2$ID,
                      ontology = "BP", method = "Wang",
                      orgdb = "org.Hs.eg.db",
                      xlab = "GO List 1", ylab = "GO List 2",
                      main = "", cex = 10)

# Comparing clusters of GO-terms ----
graph <- build_hierarchical_graph(go_list1 = data2[data2$Category == "BP",]$ID,
                                  nb_lists = "single",
                                  go_sim_object = NULL)
expanded_graph <- expand_graph(graph = graph,
                               go_list1 = data2[data2$Category == "BP",]$ID,
                               nb_lists = "single")
final_graph <- retain_ancestors_above_input_terms(graph = expanded_graph,
                                                  go_list1 = data2[data2$Category == "BP",]$ID,
                                                  nb_lists = "single")

# Semantic similarity clustering ----
similarity_matrix <- calculate_wang(graph = final_graph,
                                    ontology = "BP",
                                    orgdb = "org.Hs.eg.db")
cluster2 <- cluster_go_terms(method_type = "similarity", method = "wang",
                            orgdb = "org.Hs.eg.db", ontology = "BP",
                            similarity_matrix = similarity_matrix,
                            graph = final_graph, nb_clusters = NULL,
                            k_range = 2:25)
# "Optimal k using Silhouette Method: 24"

compare_clusters(cluster_list1 = cluster, cluster_list2 = cluster2,
                 ontology = "BP", OrgDb = "org.Hs.eg.db",
                 method = "Wang", combine = "BMA", plot = T,
                 low = "white", high = "red3",
                 labs = c("List 1", "List 2"), cex = 3)
