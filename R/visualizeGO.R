#' Visualize Hierarchy GO-terms Relationship
#'
#' This is the final function, that computes all the information in a single plot with the hierarchical relationship between the input GO-terms.
#'
#' @details
#' This function has a lot of arguments to make the plot highly customizable. I recommend ‘playing’ with the parameters to see which plot best suits what you are looking for 😊
#'
#' @param go_list1 List of GO-terms of interest.
#' @param go_list2 If enabled, second list of GO-terms to compare with the first one.
#' @param nb_lists Allows to select whenever you are gonna use one ("single") or two ("double") lists.
#' @param go_sim_object GO-terms similarity relationships saved as environment.
#' @param shape1 Shape of the nodes in the plot. By default is "circle".
#' @param shape2 If two lists are used, the shape of the nodes for the second list.
#' @param ontology Gene Ontology category to use (could be "BP" for Biological Process,
#' "CC" for Cellular Component or "MF" for "Molecular Function").
#' @param simplification If you want to filter out the GO-terms that are going to be in the plot. Default is FALSE.
#' @param min_node_size Minimum size of the nodes less connected with others. If is NULL, the size is calculated.
#' @param max_node_size Maximum size of the nodes more connected with the others. If is NULL, the size is calculated.
#' @param layout Arrangement of the nodes in the network. Can be ‘tree’, ‘kk’ or ‘fr’
#' if the chosen arrangement is tree-like, or using the Kamada-Kawai or Fruchterman-Reingold
#' algorithms respectively. For more information see the bullet points in the R \href{https://igraph.org/r/}{igraph} package.
#' @param clustering Argument indicating whether nodes are to be clustered. Default is TRUE.
#' @param clusters If the above argument is set to ‘TRUE’, result of the \code{\link[visualizeGO]{cluster_go_terms}}
#' function where the clustering information is stored.
#' @param col_palette Colors to use for each cluster.
#' @param verbose You choose the number of messages you want to be displayed on the console.
#' Different from the rest of the packages; choose between ‘none’, so that no message is produced,
#' ‘some’ if you want to generate messages about how the process is going, or ‘all’ if you want that,
#' in addition to the messages that indicate how the whole process is going, you also get intermediate tables with information.
#' @param legend Whether you want to display the legend. Default = TRUE
#' @param labs Name of each of the GO-terms lists
#' @param save_plot Set this option to "TRUE" if you want to save the plot as a PNG file. Default is FALSE.
#' @param PNG If the "save_plot" option is set to "TRUE", name of the PNG file generated.
#'
#' @return This function returns the final graph, with the colours of the nodes depending
#' on the clustering and the labels as numbers corresponding to the GO IDs in the data frame;
#' and also a table with the correspondence of the number appearing in the plot and the GO ID.
#' @export

visualizeGO <- function(cluster,
                        selected_cluster = NULL,
                        ontology = c("BP", "CC", "MF"),
                        shape1 = "circle",
                        shape2 = NULL,
                        min_node_size = NULL,
                        max_node_size = NULL,
                        layout = c("tree", "kk", "fr"),
                        col_palette = NULL,
                        verbose = c("all", "none", "some"),
                        legend = T, save_plot = F, PNG = NULL){

  # Step 0: Manage the verbose ----
  if (!verbose %in% c("all", "none", "some")) {
    stop("Verbose must be either 'all', 'some' or 'none'.")
  }

  # Step 1: Validate inputs ----
  if (verbose != "none") cat("Validating inputs...\n")
  if (is.null(cluster)) {
    stop("cluster cannot be NULL or empty.")
  }

  if (is.null(cluster) || !"clusters" %in% names(cluster) || !"graph" %in% names(cluster)) {
    stop("Invalid cluster result provided. It should contain 'clusters' and 'graph' components.")
  }

  cluster_assignments <- cluster$clusters
  graph <- cluster$graph

  # Step 2: Validate vertex names ----
  if (verbose != "none") cat("Validating vertex names...\n")
  if (is.null(V(graph)$name) || any(is.na(V(graph)$name)) || any(duplicated(V(graph)$name))) {
    stop("Invalid vertex names detected. Ensure all vertices have unique, non-NA names.")
  }

  # Step 3: Filter nodes by selected cluster(s) if enabled ----
  if (!is.null(selected_cluster)) {
    # Keep only the nodes belonging to the selected cluster(s)
    selected_nodes <- V(graph)[V(graph)$cluster %in% selected_cluster]

    # Subset the graph to only include the selected nodes
    graph <- induced_subgraph(graph, selected_nodes)
  }

  # Step 4: Determine cluster colors ----
  unique_clusters <- unique(V(graph)$cluster)
  unique_clusters <- unique_clusters[!is.na(unique_clusters)]  # Exclude NA values

  # Ensure col_palette matches the number of unique clusters
  if (length(col_palette) != length(unique_clusters)) {
    stop("The length of 'col_palette' must match the number of unique clusters.")
  }

  # First, ensure the length of col_palette matches the number of selected nodes
  if (length(col_palette) == 1) {
    # If only one color is provided, apply it to all selected nodes
    V(graph)$color <- col_palette[1]
    } else if (is.null(col_palette) || length(col_palette) == 0) {
      cat("Generating a dynamic color palette.\n")
      col_palette <- generate_pastel_colors(length(unique(V(graph)$cluster)))
    }

  cluster_colors <- setNames(col_palette, unique_clusters)

  V(graph)$color <- ifelse(
    is.na(V(graph)$cluster),
    "white",  # Default color for nodes without a cluster
    cluster_colors[as.character(V(graph)$cluster)])

  # Step 5: Assign shapes to nodes ----
  if (verbose != "none") cat("Assign shapes to nodes ...\n")
  nb_lists <- dplyr::case_when(
    "List2" %in% unique(V(graph)$list) ~ "double",
    TRUE ~ "single")

  if (nb_lists == "single"){
    V(graph)$shape <- ifelse(!is.na(V(graph)$list), shape1, "circle")
  } else {
    if (is.null(shape2)){
      stop("You must chose the shape for the second GO terms list ('shape2' argument).")
    } else {
      V(graph)$shape <- dplyr::case_when(
        nb_lists == "single" ~ ifelse(!is.na(V(graph)$list), shape1, "circle"),
        nb_lists == "double" ~ ifelse(V(graph)$list == "List2", shape2,
                                      !is.na(V(graph)$list), shape1, "circle"))  # Default shape for other nodes
    }
  }

  if(verbose == "all"){
    cat("Vertex Names and Shapes:\n")
    print(data.frame(Name = V(graph)$name, Shape = V(graph)$shape))
  }

  # Step 6: Scale node sizes ----
  if (verbose != "none") cat("Scale nodes sizes ...\n")
  node_sizes <- degree(graph)
  # If no user-defined size range, calculate it from the node degrees
  if (is.null(min_node_size)) {
    min_node_size <- min(node_sizes)
    cat("Minimum Node Size: ", min_node_size, "\n")
  }
  if (is.null(max_node_size)) {
    max_node_size <- max(node_sizes)
    cat("Maximum Node Size: ", max_node_size, "\n")
  }

  # Scale node sizes based on the min and max values defined by the user
  scaled_node_sizes <- (node_sizes -
                          min(node_sizes)) /
    (max(node_sizes) - min(node_sizes)) *
    (max_node_size - min_node_size) +
    min_node_size
  V(graph)$size <- scaled_node_sizes

  # Step 7: Layout and transformation ----
  if (verbose != "none") cat("Transforming layout ...\n")
  layout <- switch(layout,
                   tree = layout_as_tree(graph, root = V(graph)[degree(graph, mode = "in") == 0]),
                   fr = layout_with_fr(graph, niter = 500, grid = "nogrid"),
                   layout_with_kk(graph, niter = 1000, kkconst = 0.2, maxiter = 15000))

  if (any(!is.finite(layout))) {
    stop("Error: The layout contains non-finite values.")
  }

  # Step 8: Transform GO IDs to numbers ----
  if (verbose != "none") cat("Transform GO IDs to numbers ...\n")
  transformation_result <- transform_go_ids_to_numbers(graph)
  graph <- transformation_result$graph
  go_id_to_description <- transformation_result$go_id_to_description

  GODescriptions <- transformation_result$GODescriptions
  # This stores the GODescriptions data frame

  # Debugging: Print go_terms and go_descriptions to check the output
  if (verbose == "all"){
    print(paste("GO Terms: ", paste(go_terms, collapse = ", ")))
    print(paste("GO Descriptions: ", paste(go_descriptions, collapse = ", ")))

    # Debugging: Print the GODescriptions data frame
    print("GODescriptions Data Frame:")
    print(GODescriptions)
  }

  # Step 9: Plot the graph ----
  if (verbose != "none") cat("Displaying graph...\n")

  if (vcount(graph) == 0 || ecount(graph) == 0) {
    stop("Error: The graph is empty. No nodes or edges to plot.")
  }

  par(mar = c(1.5, 1.5, 1.5, 1.5))

  plot(
    graph,
    layout = layout,
    vertex.frame.color = "black",
    vertex.label = V(graph)$name,
    vertex.size = V(graph)$size,
    vertex.label.color = "black",
    vertex.color = V(graph)$color,
    vertex.frame.color = "black",
    vertex.label.cex = 0.7,
    vertex.label.family = "sans",  # Set the font family to "sans"
    vertex.shape = V(graph)$shape,
    edge.arrow.size = 0.3,
    edge.color = "darkgray",
    main = paste("Hierarchical GO:", ontology, " Graph", sep = ""),
    rescale = TRUE, # Allow the graph to scale to fit the available space
    margin = 0          # Remove additional margins from the plot
  )

  # Step 10: Add legend ----
  if(isTRUE(legend)){
    # Add a legend for clusters
    if (!is.null(V(graph)$cluster)) {
      cluster_labels <- paste("Cluster", unique(V(graph)$cluster))  # Create cluster labels
      cluster_labels <- ifelse(cluster_labels == "Cluster NA", "No Cluster", cluster_labels)
      cluster_legend_colors <- cluster_colors[unique(as.character(V(graph)$cluster))]  # Match colors
      legend("topleft",
             legend = cluster_labels,
             fill = cluster_legend_colors,
             bty = "n", title.font = 2,
             cex = 0.8, title = "Clusters",
             inset = c(0.02, 0.001))
    }

    # Add a legend for node shapes
    # Map shapes to pch values
    shape_to_pch <- ifelse(
      V(graph)$shape == "circle", 21,  # `pch = 21` corresponds to filled circle
      ifelse(V(graph)$shape == "square", 22,  # `pch = 22` corresponds to filled square
             ifelse(V(graph)$shape == "rectangle", 23,  # `pch = 23` corresponds to filled rectangle
                    NA)))  # Default to NA for unsupported shapes

    # Use preprocessed pch in the legend
    unique_shapes <- unique(V(graph)$shape)
    legend_shapes <- unique(shape_to_pch[!is.na(shape_to_pch)])
    if(length(unique(V(graph)$cluster)) <= 15){
      if(nb_lists == "single"){
        legend("topleft",
               legend = c(labs, "Other Terms"),
               pch = legend_shapes, # Extract unique pch values for the legend
               bty = "n", title.font = 2, cex = 0.8,
               title = "Node Origin", xjust = 1, inset = c(0.02, 0.7))
      } else {
        legend("topleft",
               legend = c(labs, "Other Terms"),
               pch = legend_shapes, # Extract unique pch values for the legend
               bty = "n", title.font = 2, cex = 0.8,
               title = "Node Origin", xjust = 1, inset = c(0.02, 0.7))
      }
    } else {
      if(nb_lists == "single"){
        legend("topright",
               legend = c(labs, "Other Terms"),
               pch = legend_shapes, # Extract unique pch values for the legend
               bty = "n", title.font = 2, cex = 0.8,
               title = "Node Origin", xjust = 1, inset = c(0.02, 0.7))
      } else {
        legend("topright",
               legend = c(labs, "Other Terms"),
               pch = legend_shapes, # Extract unique pch values for the legend
               bty = "n", title.font = 2, cex = 0.8,
               title = "Node Origin", xjust = 1, inset = c(0.02, 0.7))
      }
    }

    # Add a legend for the node sizes (degree of connectivity)
    if(length(unique(V(graph)$cluster)) <= 15){
      legend("topleft", legend = c("Low Connectivity",
                                   "High Connectivity"),
             pch = 21, pt.bg = "lightgray", title = "Degree of connectivity",
             pt.cex = c(min(scaled_node_sizes),
                        max(scaled_node_sizes)/3),
             bty = "n", cex = 0.8, title.font = 2,
             xjust = 1, inset = c(0.01, 0.85))
    } else {
      legend("topright", legend = c("Low Connectivity",
                                    "High Connectivity"),
             pch = 21, pt.bg = "lightgray", title = "Degree of connectivity",
             pt.cex = c(min(scaled_node_sizes),
                        max(scaled_node_sizes)/3),
             bty = "n", cex = 0.8, title.font = 2,
             xjust = 1, inset = c(-0.035, 0.85))
    }

  }

  # Step 11: Save plot (if enabled) ----
  if(isTRUE(save_plot)) {
    if (verbose != "none") cat("Saving plot as", PNG, "...\n")

    # Ensure the file extension is included in the filename (if not already)
    if (!grepl("\\.png$", PNG)) {
      file_name <- paste0(PNG, ".png")  # Default to .png if no extension is provided
    }

    # Open a PNG device to save the plot with high resolution
    png(PNG, width = 1427, height = 674, res = 100)

    # PLOT ----
    par(mar = c(1.5, 1.5, 1.5, 1.5))
    plot(
      graph,
      layout = layout,
      vertex.frame.color = "black",
      vertex.label = V(graph)$name,
      vertex.size = V(graph)$size,
      vertex.label.color = "black",
      vertex.color = V(graph)$color,
      vertex.frame.color = "black",
      vertex.label.cex = 0.7,
      vertex.label.family = "sans",  # Set the font family to "sans"
      vertex.shape = V(graph)$shape,
      edge.arrow.size = 0.5,
      edge.color = "darkgray",
      main = paste("Hierarchical GO:", ontology, " Graph", sep = ""),
      rescale = TRUE, # Allow the graph to scale to fit the available space
      margin = 0          # Remove additional margins from the plot
    )
    # LEGEND ----
    if(isTRUE(legend)){
      # Add a legend for clusters
      legend("topleft",
             legend = cluster_labels,
             fill = cluster_legend_colors,
             bty = "n", title.font = 2,
             cex = 0.8, title = "Clusters",
             inset = c(0.02, 0.001))

      # Add legend for node origin (shapes)
      if(nb_lists == "single"){
        legend("topleft",
               legend = c(labs, "Other Terms"),
               pch = legend_shapes, # Extract unique pch values for the legend
               bty = "n", title.font = 2, cex = 0.8,
               title = "Node Origin", xjust = 1, inset = c(0.02, 0.7))
      } else {
        legend("topleft",
               legend = c(labs, "Other Terms"),
               pch = legend_shapes, # Extract unique pch values for the legend
               bty = "n", title.font = 2, cex = 0.8,
               title = "Node Origin", xjust = 1, inset = c(0.02, 0.7))
      }

      # Add a legend for the node sizes (degree of connectivity)
      legend("topleft", legend = c("Low Connectivity",
                                   "High Connectivity"),
             pch = 21, pt.bg = "lightgray", title = "Degree of connectivity",
             pt.cex = c(min(scaled_node_sizes),
                        max(scaled_node_sizes)/3),
             bty = "n", cex = 0.8, title.font = 2,
             xjust = 1, inset = c(0.022, 0.8))
    }

    dev.off()
    if (verbose != "none") cat("Plot saved successfully as", PNG, "\n")
  }
  # Step 12: Return the GODescriptions ----
  if(verbose == "all"){
    return(GODescriptions)
  }
}
