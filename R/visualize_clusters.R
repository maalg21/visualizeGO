#' Visualize Clusters of GO Terms
#'
#' The visualize_clusters function is designed to prepare a graph for visualization
#' by assigning unique colors to clusters in the graph based on their memberships.
#' This function is primarily used to prepare a graph for visualization,
#' ensuring that clusters are visually distinguishable based on their assigned colors.
#'
#' @param graph A graph object that connects the input GO-terms with their parents and children terms groupped in clusters.
#' @param col_palette If no color palette is provided (NULL), the function generates a default palette using rainbow, with one color for each unique cluster in the graph.
#' @return Each node in the graph is assigned a color based on its cluster membership. The cluster colors are mapped from the color palette. The function returns a named list, where each cluster is mapped to its assigned color.
#' @export

visualize_clusters <- function(graph, col_palette = NULL) {

  # Default color palette for clusters
  if (is.null(col_palette)) {
    col_palette <- rainbow(length(unique(V(graph)$cluster)))
  }

  # Assign colors to clusters
  cluster_colors <- col_palette[V(graph)$cluster]

  # Set node colors in the graph
  V(graph)$color <- cluster_colors

  # Return a named list of colors for each cluster
  cluster_color_list <- setNames(col_palette, unique(V(graph)$cluster))

  return(cluster_color_list)  # Return the list of colors
}
