devtools::install_github("maalg21/visualizeGO", force = T)
library(visualizeGO)

data <- as.data.frame(readxl::read_xlsx(system.file("extdata", "GOTerms.xlsx",
                                                    package = "visualizeGO")))
data <- data[,c("Category", "ID")]
GO_BP <- data[data$Category == "BP",]

graph <- build_hierarchical_graph(go_list1 = GO_BP$ID, nb_lists = "single",
                                  go_sim_object = NULL)
expanded_graph <- expand_graph(graph = graph,
                               go_list1 = GO_BP$ID,
                               nb_lists = "single")
final_graph <- retain_ancestors_above_input_terms(graph = expanded_graph,
                                                  go_list1 = GO_BP$ID,
                                                  nb_lists = "single")

similarity_matrix <- calculate_wang(graph = final_graph,
                                    ontology = "BP",
                                    orgdb = "org.Hs.eg.db")
cluster <- cluster_go_terms(method_type = "similarity", method = "wang",
                            orgdb = "org.Hs.eg.db", ontology = "BP", similarity_matrix = similarity_matrix,
                            graph = final_graph, nb_clusters = NULL, k_range = 2:10)

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
                       verbose = "some", # Just to know more about the process that is ocurring
                       save_plot = T, PNG = "plot1.png")

filter_and_visualize_cluster <- function(cluster_go_result = cluster,
                                         selected_cluster = 5,
                                         go_sim_object = NULL,
                                         shape = "circle",
                                         ontology = "BP",
                                         layout = "tree",
                                         col_palette = "#D9C9FF",
                                         min_node_size = 1,
                                         max_node_size = 10,
                                         save_plot = FALSE,
                                         PNG = NULL,
                                         verbose = "none")
