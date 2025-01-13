#' Visualize Clusters GO-terms Relationship
#'
#' This function  computes all the information in a single plot with the hierarchical relationship between the input GO-terms from a single selected cluster.
#'
#' @param clusters Result of the \code{\link[visualizeGO]{cluster_go_terms}} function where the clustering information is stored.
#' @param selected_cluster Cluster that we want to analyse more in detail.
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

filter_and_visualize_cluster <- function(clusters,
                                         selected_cluster = NULL,
                                         ontology = c("BP", "CC", "MF"),
                                         min_node_size = NULL,
                                         max_node_size = NULL,
                                         layout = c("tree", "kk", "fr"),
                                         col_palette = NULL,
                                         verbose = c("all", "none", "some"),
                                         save_plot = FALSE, PNG = NULL) {

  # Step 1: Validate inputs ----
  if (!verbose %in% c("all", "none", "some")) {
    stop("Verbose must be either 'all', 'some' or 'none'.")
  }

  if (is.null(cluster_result) || !"clusters" %in% names(cluster_result) || !"graph" %in% names(cluster_result)) {
    stop("Invalid cluster result provided. It should contain 'clusters' and 'graph' components.")
  }

  cluster_assignments <- clusters$clusters
  graph <- clusters$graph

  # Step 2: Filter nodes by selected cluster(s) ----
  if (is.null(selected_cluster)) {
    stop("Please specify the selected cluster.")
  }

  # Keep only the nodes belonging to the selected cluster(s)
  selected_nodes <- V(graph)[V(graph)$cluster %in% selected_cluster]

  # Subset the graph to only include the selected nodes
  subgraph <- induced_subgraph(graph, selected_nodes)

  # Step 3: Prepare the graph for visualization ----
  if (verbose != "none") cat("Building the filtered graph...\n")

  # Assign node colors based on clusters
  if (is.null(col_palette) || length(col_palette) == 0) {
    col_palette <- generate_pastel_colors(length(unique(V(subgraph)$cluster)))
  }

  cluster_colors <- setNames(col_palette, unique(V(subgraph)$cluster))
  V(subgraph)$color <- cluster_colors[as.character(V(subgraph)$cluster)]

  # Assign shapes to nodes (you can customize shapes here)
  V(subgraph)$shape <- "circle"

  # Scale node sizes (based on degree or any other criteria)
  node_sizes <- degree(subgraph)

  # If no user-defined size range, calculate it from the node degrees
  if (is.null(min_node_size)) {
    min_node_size <- min(node_sizes)
  }
  if (is.null(max_node_size)) {
    max_node_size <- max(node_sizes)
  }

  # Scale node sizes based on the min and max values
  scaled_node_sizes <- (node_sizes - min(node_sizes)) / (max(node_sizes) - min(node_sizes)) *
    (max_node_size - min_node_size) + min_node_size
  V(subgraph)$size <- scaled_node_sizes

  # Step 4: Layout transformation ----
  layout <- match.arg(layout)
  layout_fun <- switch(layout,
                       tree = layout_as_tree(subgraph, root = V(subgraph)[degree(subgraph, mode = "in") == 0]),
                       fr = layout_with_fr(subgraph, niter = 500, grid = "nogrid"),
                       kk = layout_with_kk(subgraph, niter = 1000, kkconst = 0.2, maxiter = 15000))

  # Step 5: GO ID and Term Relationship (Create a data frame) ----
  # Obtain the descriptions of the GO IDs
  go_ids <- V(subgraph)$name
  go_descriptions <- AnnotationDbi::Term(go_ids)

  go_id_to_description <- data.frame(
    GO_ID = go_ids,
    Description = go_descriptions
  )

  # Print the relationship data frame (optional, based on verbose)
  if (verbose == "all") {
    print(go_id_to_description)
  }

  # Step 6: Plot the graph ----
  if (verbose != "none") cat("Displaying the graph...\n")
  plot(subgraph,
       layout = layout_fun,
       vertex.frame.color = "black",
       vertex.label = V(subgraph)$name,
       vertex.size = V(subgraph)$size,
       vertex.label.color = "black",
       vertex.color = V(subgraph)$color,
       vertex.frame.color = "black",
       vertex.label.cex = 0.7,
       vertex.label.family = "sans",  # Set the font family to "sans"
       vertex.shape = V(subgraph)$shape,
       edge.arrow.size = 0.5,
       edge.color = "darkgray",
       main = paste("Cluster Visualization: ", paste(selected_cluster, collapse = ", "), " Ontology: ", ontology, sep = ""),
       rescale = TRUE)

  # Step 7: Add legend ----
  if (!is.null(V(subgraph)$cluster)) {
    cluster_labels <- paste("Cluster", unique(V(subgraph)$cluster))
    cluster_labels <- ifelse(cluster_labels == "Cluster NA", "No Cluster", cluster_labels)
    cluster_legend_colors <- cluster_colors[unique(as.character(V(subgraph)$cluster))]
    legend("topleft",
           legend = cluster_labels,
           fill = cluster_legend_colors,
           bty = "n", title.font = 2,
           cex = 0.8, title = "Clusters",
           inset = c(0.001, 0.05))
  }

  # Add GO IDs and Terms below the clusters in the legend
  go_terms_table <- go_id_to_description[go_id_to_description$GO_ID %in% V(subgraph)$name, ]
  go_terms_text <- apply(go_terms_table, 1, function(x) paste(x[1], ":", x[2]))
  legend("topleft",
         legend = go_terms_text,
         bty = "n", inset = c(0.001, -0.1), cex = 0.6, title = "GO IDs and Terms")

  legend("topright", legend = c("Low Connectivity",
                                "High Connectivity"),
         pch = 21, pt.bg = "lightgray", title = "Degree of connectivity",
         pt.cex = c(1, max(scaled_node_sizes)/10),
         bty = "n", cex = 0.8, title.font = 2,
         xjust = 1, inset = c(0.00009, 0.8))

  # Step 8: Save plot (if enabled) ----
  if (isTRUE(save_plot)) {
    if (verbose != "none") cat("Saving plot as", PNG, "...\n")

    if (!grepl("\\.png$", PNG)) {
      PNG <- paste0(PNG, ".png")  # Default to .png if no extension is provided
    }

    # Set high resolution for the plot (e.g., 300 DPI)
    dpi <- 300

    # Open a PNG device to save the plot with high resolution
    png(PNG, width = 4000, height = 3000, res = dpi)

    # PLOT ----
    plot(subgraph,
         layout = layout_fun,
         vertex.frame.color = "black",
         vertex.label = V(subgraph)$name,
         vertex.size = V(subgraph)$size,
         vertex.label.color = "black",
         vertex.color = V(subgraph)$color,
         vertex.frame.color = "black",
         vertex.label.cex = 0.7,
         vertex.label.family = "sans",
         vertex.shape = V(subgraph)$shape,
         edge.arrow.size = 0.5,
         edge.color = "darkgray",
         main = paste("Cluster ", paste(selected_cluster, collapse = ", "), "\nGene Ontology ", ontology, sep = ""),
         rescale = TRUE)

    # LEGEND ----
    if (!is.null(V(subgraph)$cluster)) {
      cluster_labels <- paste("Cluster", unique(V(subgraph)$cluster))
      cluster_labels <- ifelse(cluster_labels == "Cluster NA", "No Cluster", cluster_labels)
      cluster_legend_colors <- cluster_colors[unique(as.character(V(subgraph)$cluster))]
      legend("topleft",
             legend = cluster_labels,
             fill = cluster_legend_colors,
             bty = "n", title.font = 2,
             cex = 0.8, title = "Clusters",
             inset = c(0.001, 0.05))
    }

    # Add GO IDs and Terms below the clusters in the legend
    go_terms_table <- go_id_to_description[go_id_to_description$GO_ID %in% V(subgraph)$name, ]
    go_terms_text <- apply(go_terms_table, 1, function(x) paste(x[1], ":", x[2]))
    legend("topleft",
           legend = go_terms_text,
           bty = "n", inset = c(0.001, -0.1), cex = 0.6, title = "GO IDs and Terms")

    legend("topright", legend = c("Low Connectivity",
                                  "High Connectivity"),
           pch = 21, pt.bg = "lightgray", title = "Degree of connectivity",
           pt.cex = c(1, max(scaled_node_sizes)/10),
           bty = "n", cex = 0.8, title.font = 2,
           xjust = 1, inset = c(0.00009, 0.8))

    dev.off()
    if (verbose != "none") cat("Plot saved successfully as", PNG, "\n")
  }
}
