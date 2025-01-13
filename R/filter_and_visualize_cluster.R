#' Visualize Clusters GO-terms Relationship
#'
#' This function  computes all the information in a single plot with the hierarchical relationship between the input GO-terms from a single selected cluster.
#'
#' @param clusters Result of the \code{\link[visualizeGO]{cluster_go_terms}} function where the clustering information is stored.
#' @param selected_cluster Cluster that we want to analyse more in detail.
#' @param go_sim_object GO-terms similarity relationships saved as environment.
#' @param shape Shape of the nodes in the plot. By default is "circle".
#' @param ontology Gene Ontology category to use (could be "BP" for Biological Process,
#' "CC" for Cellular Component or "MF" for "Molecular Function").
#' @param min_node_size Minimum size of the nodes less connected with others. If is NULL, the size is calculated.
#' @param max_node_size Maximum size of the nodes more connected with the others. If is NULL, the size is calculated.
#' @param layout Arrangement of the nodes in the network. Can be ‘tree’, ‘kk’ or ‘fr’
#' if the chosen arrangement is tree-like, or using the Kamada-Kawai or Fruchterman-Reingold
#' algorithms respectively. For more information see the bullet points in the R \href{https://igraph.org/r/}{igraph} package.
#' @param col_palette Colors to use for the selected cluster.
#' @param verbose You choose the number of messages you want to be displayed on the console.
#' Different from the rest of the packages; choose between ‘none’, so that no message is produced,
#' ‘some’ if you want to generate messages about how the process is going, or ‘all’ if you want that,
#' in addition to the messages that indicate how the whole process is going, you also get intermediate tables with information.
#' @param save_plot Set this option to "TRUE" if you want to save the plot as a PNG file. Default is FALSE.
#' @param PNG If the "save_plot" option is set to "TRUE", name of the PNG file generated.
#'
#' @return This function returns the graph of the selected cluster.
#' @export

filter_and_visualize_cluster <- function(cluster_go_result,
                                         selected_cluster,
                                         go_sim_object = NULL,
                                         shape = "circle",
                                         ontology = c("BP", "CC", "MF"),
                                         layout = c("tree", "kk", "fr"),
                                         col_palette, min_node_size,
                                         max_node_size,
                                         save_plot = FALSE,
                                         PNG = NULL,
                                         verbose = c("all", "none", "some")) {

  # Validate verbose input
  verbose <- match.arg(verbose)

  if (verbose != "none") cat("Filtering for selected cluster...\n")

  # Validate that the result of cluster_go_terms contains ‘graph’ and ‘clusters’.
  if (!"graph" %in% names(cluster_go_result) || !"clusters" %in% names(cluster_go_result)) {
    stop("The cluster_go_result must contain both 'graph' and 'clusters'.")
  }

  graph <- cluster_go_result$graph
  clusters <- cluster_go_result$clusters

  # Filter out nodes belonging to the selected cluster
  selected_nodes <- names(clusters[clusters == selected_cluster])

  # Create a subnetwork with only the nodes of the selected cluster
  subgraph <- induced_subgraph(graph, vids = selected_nodes)

  # Display the subnetwork using visualize_go_hierarchy
  if (verbose != "none") cat("Visualizing cluster...\n")

  visualize_go_hierarchy(go_list1 = selected_nodes,
                         go_sim_object = go_sim_object,
                         nb_lists = "single",
                         shape1 = shape,
                         ontology = ontology,
                         layout = layout,
                         clustering = FALSE,
                         col_palette = col_palette,
                         min_node_size = min_node_size,
                         max_node_size = max_node_size,
                         save_plot = save_plot,
                         PNG = PNG,
                         verbose = verbose)

  # Return the sub-graph for further analysis if necessary.
  return(subgraph)
}
