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
#' @param legend Whether you want to display de legend. Default = TRUE
#' @param labs Name of each of the GO-terms lists
#' @param save_plot Set this option to "TRUE" if you want to save the plot as a PNG file. Default is FALSE.
#' @param PNG If the "save_plot" option is set to "TRUE", name of the PNG file generated.
#'
#' @return This function returns the final graph, with the colours of the nodes depending
#' on the clustering and the labels as numbers corresponding to the GO IDs in the data frame;
#' and also a table with the correspondence of the number appearing in the plot and the GO ID.
#' @export
visualize_go_hierarchy <- function(go_list1, go_list2 = NULL,
                                   nb_lists = c("single", "double"),
                                   go_sim_object = NULL,
                                   shape1 = "circle",
                                   shape2 = NULL,
                                   ontology = c("BP", "CC", "MF"),
                                   simplification = F,
                                   min_node_size = NULL,
                                   max_node_size = NULL,
                                   layout = c("tree", "kk", "fr"),
                                   clustering = T,
                                   clusters = NULL,
                                   col_palette = NULL,
                                   verbose = c("all", "none", "some"),
                                   legend = T, labs = NULL,
                                   save_plot = F,
                                   PNG = NULL
){

  # Step 0: Manage the verbose ----
  if (!verbose %in% c("all", "none", "some")) {
    stop("Verbose must be either 'all', 'some' or 'none'.")
  }
  # Step 1: Validate inputs ----
  if (verbose != "none") cat("Validating inputs...\n")

  if(is.null(go_sim_object)){
    # Get go_sim_object
    file_name <- c("go_term_database_all_ontologies.RData")
    file <- system.file("extdata", file_name, package = "visualizeGO")

    # Check if the file exists in the specified directory
    if (!file.exists(file)) {
      stop(paste("Error: The file", file_name, "does not exist in the directory or it's corrupted."))
    } else {
      message("File found: ", file_name)
      load(file)
      go_sim_object <- GOSimEnv
    }
  }

  if (is.null(go_list1) || length(go_list1) == 0) {
    stop("go_list1 cannot be NULL or empty.")
  }

  # Step 2: Build the graph ----
  if (verbose != "none") cat("Building GO graph...\n")
  if(nb_lists == "single"){
    # Build the graph
    graph <- build_hierarchical_graph(go_list1)
    if(verbose != "none") cat(paste(vcount(graph), " vertices considered\n", sep = ""))
  } else if(nb_lists == "double") {
    # Build the graph
    graph <- build_hierarchical_graph(go_list1, go_list2)
    if(verbose != "none") cat(paste(vcount(graph), " vertices considered\n", sep = ""))
  }

  # Step 3: Validate vertex names ----
  if (verbose != "none") cat("Validating vertex names...\n")
  if (is.null(V(graph)$name) || any(is.na(V(graph)$name)) || any(duplicated(V(graph)$name))) {
    stop("Invalid vertex names detected. Ensure all vertices have unique, non-NA names.")
  }

  # Step 4: Simplification (if enabled) ----
  if (verbose != "none" && isTRUE(simplification)) cat("Simplifying GO graph...\n")
  if (isTRUE(simplification)) {
    if (nb_lists == "single"){
      graph <- expand_graph(graph, go_list1, nb_lists)
      graph <- retain_ancestors_above_input_terms(graph, go_list1, nb_lists)
    } else if (nb_lists == "double") {
      graph <- expand_graph(graph, go_list1, go_list2, nb_lists)
      graph <- retain_ancestors_above_input_terms(graph, go_list1, go_list2, nb_lists)
    }
    if(verbose != "none") cat(paste(vcount(graph), " vertices after the filtering\n", sep = ""))
  } else {
    cat("You have not filtered the connections in your graph.\nIt is possible that, if you have many nodes, the graph is not explanatory.\n")
  }

  # Re-check vertex names after simplification
  if (is.null(V(graph)$name) || any(is.na(V(graph)$name)) || any(duplicated(V(graph)$name))) {
    stop("Invalid vertex names detected after simplification. Ensure unique, non-NA names.")
  }

  # Step 4: Clustering (if enabled) ----
  # If clusters is NULL, display a message or perform an alternative action
  if (verbose != "none" && is.null(clusters)) {
    message("No clusters were provided, generating the graph without clustering.")

    # Assign colors to all nodes based on col_palette
    if (is.null(col_palette) || length(col_palette) == 0) {
      # If no col_palette is provided, generate a default color
      message("No color palette provided. Using default color palette.")
      col_palette <- generate_pastel_colors(1)  # Single color for all nodes
    }

    # Ensure the length of col_palette matches the number of selected nodes
    if (length(col_palette) == 1) {
      # If only one color is provided, apply it to all selected nodes
      V(graph)$color <- col_palette[1]
    } else if (length(col_palette) >= length(V(graph))) {
      # If col_palette has enough colors, assign them sequentially
      V(graph)$color <- col_palette[seq_along(V(graph))]
    } else {
      stop("The provided 'col_palette' has insufficient colors for the selected nodes.")
    }

  } else {
    message("Clusters were provided, generating the graph with groupings.")
  }

  if (verbose != "none" && !is.null(clusters)) cat("Applying clustering...\n")
  if (!is.null(clusters) && isTRUE(clustering)) { # Cluster inclusion logic

    # Extract the actual cluster assignments
    if (!"clusters" %in% names(clusters) || !"graph" %in% names(clusters)) {
      stop("The 'clusters' input must contain both 'graph' and 'clusters' components.")
    }

    cluster_assignments <- clusters$clusters
    graph <- clusters$graph

    # Validate that clusters match graph vertices
    if (!all(V(graph)$name %in% names(cluster_assignments))) {
      unmatched_vertices <- setdiff(V(graph)$name, names(cluster_assignments))
      stop(paste("Cluster information does not match the graph vertices.\nThe following vertices are missing in cluster assignments:",
                 paste(unmatched_vertices, collapse = ", ")))
    }

    # Assign clusters to graph vertices
    V(graph)$cluster <- cluster_assignments[V(graph)$name]

    # Determine cluster colors
    unique_clusters <- unique(V(graph)$cluster)
    unique_clusters <- unique_clusters[!is.na(unique_clusters)]  # Exclude NA values

    # Handle NA clusters explicitly
    if (any(is.na(V(graph)$cluster))) {
      cat("NA clusters detected. Nodes without cluster assignments will be displayed with 'white' color.\n")
    }

    # Generate colors dynamically if col_palette is NULL or insufficient
    if (is.null(col_palette) || length(col_palette) < length(unique_clusters)) {
      cat("Generating a dynamic color palette.\n")
      col_palette <- generate_pastel_colors(length(unique_clusters))
    }

    # Ensure col_palette matches the number of unique clusters
    if (length(col_palette) != length(unique_clusters)) {
      stop("The length of 'col_palette' must match the number of unique clusters.")
    }

    cluster_colors <- setNames(col_palette, unique_clusters)

    V(graph)$color <- ifelse(
      is.na(V(graph)$cluster),
      "white",  # Default color for nodes without a cluster
      cluster_colors[as.character(V(graph)$cluster)])
  }

  # Step 5: Assign shapes to nodes ----
  if (verbose != "none") cat("Assign shapes to nodes ...\n")
  if (nb_lists == "single"){
    V(graph)$shape <- shape1  # All nodes have the same shape for a single list
  } else {
    if (is.null(shape2)){
      stop("You must chose the shape for the second GO terms list ('shape2' argument).")
    } else {
      V(graph)$shape <- ifelse(
        V(graph)$name %in% go_list1, shape1,
        ifelse(V(graph)$name %in% go_list2, shape2, "circle"))  # Default shape for other nodes
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
    rescale = TRUE # Allow the graph to scale to fit the available space
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
             inset = c(0.001, 0.05))
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
    if(nb_lists == "single"){
      legend("topright",
             legend = c(labs, "Other Terms"),
             pch = legend_shapes, # Extract unique pch values for the legend
             bty = "n", title.font = 2, cex = 0.8,
             title = "Node Origin", xjust = 1, inset = c(0.035, 0.6))
    } else {
      legend("topright",
             legend = c(labs, "Other Terms"),
             pch = legend_shapes, # Extract unique pch values for the legend
             bty = "n", title.font = 2, cex = 0.8,
             title = "Node Origin", xjust = 1, inset = c(0.035, 0.6))
    }

    # Add a legend for the node sizes (degree of connectivity)
    legend("topright", legend = c("Low Connectivity",
                                  "High Connectivity"),
           pch = 21, pt.bg = "lightgray", title = "Degree of connectivity",
           pt.cex = c(min_node_size, max_node_size),
           bty = "n", cex = 0.8, title.font = 2,
           xjust = 1, inset = c(0.00009, 0.8))
  }

  # Step 11: Save plot (if enabled) ----
  if(isTRUE(save_plot)) {
    if (verbose != "none") cat("Saving plot as", PNG, "...\n")

    # Ensure the file extension is included in the filename (if not already)
    if (!grepl("\\.png$", PNG)) {
      file_name <- paste0(PNG, ".png")  # Default to .png if no extension is provided
    }

    # Set high resolution for the plot (e.g., 300 DPI)
    dpi <- 300  # High resolution

    # Open a PNG device to save the plot with high resolution
    png(PNG, width = 4000, height = 3000, res = dpi)

    # PLOT ----
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
      rescale = TRUE # Allow the graph to scale to fit the available space
    )
    # LEGEND ----
    if(isTRUE(legend)){
      legend("topleft",
             legend = cluster_labels,
             fill = cluster_legend_colors,
             bty = "n", title.font = 2,
             cex = 0.8, title = "Clusters",
             inset = c(0.001, 0.05))

      if(nb_lists == "single"){
        legend("topright",
               legend = c(labs, "Other Terms"),
               pch = legend_shapes, # Extract unique pch values for the legend
               bty = "n", title.font = 2, cex = 0.8,
               title = "Node Origin", xjust = 1, inset = c(0.035, 0.6))
      } else {
        legend("topright",
               legend = c(labs, "Other Terms"),
               pch = legend_shapes, # Extract unique pch values for the legend
               bty = "n", title.font = 2, cex = 0.8,
               title = "Node Origin", xjust = 1, inset = c(0.035, 0.6))
      }

      legend("topright", legend = c("Low Connectivity",
                                    "High Connectivity"),
             pch = 21, pt.bg = "lightgray", title = "Degree of connectivity",
             pt.cex = c(min(scales_node_sizes), max(scaled_node_sizes)),
             bty = "n", cex = 0.8, title.font = 2,
             xjust = 1, inset = c(0.00009, 0.8))
    }

    dev.off()
    if (verbose != "none") cat("Plot saved successfully as", PNG, "\n")
  }
  # Step 12: Return the GODescriptions ----
  if(verbose == "all"){
    return(GODescriptions)
  }
}
