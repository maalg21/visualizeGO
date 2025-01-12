#' Generate Clusters Table
#'
#' The generate_cluster_table function creates a detailed and visually appealing table
#' summarizing the clustering results of Gene Ontology (GO) terms. The function can work
#' with both similarity-based clustering and network-based clustering outputs.
#'
#' This function is useful for reporting clustering results in a clear, interactive,
#' and visually appealing format, ideal for presentations or publications.
#'
#' @description
#' A formatted HTML table summarizing: (1) Cluster ID, (2) Representative Pathway
#' (for similarity-based clustering) and (3) GO IDs within each cluster.
#' Adds row-specific colors to highlight clusters visually.
#'
#' @param cluster_output Result from "cluster_go_terms".
#' @param method_type Method used to perform the clustering; either "similarity" or "network".
#' @param similarity_matrix If the method used was "similarity"; similarity matrix used to perform the clustering.
#' @param col_palette Colors to use for each cluster.
#' @param text_color Color of the text in the table.
#' @return Combines cluster data into a single data frame. Assigns colors to rows based on cluster memberships using the provided or default color palette.
#' @seealso [cluster_go_terms()]
#' @export

generate_cluster_table <- function(cluster_output,
                                   method_type = c("similarity", "network"),
                                   similarity_matrix = NULL,
                                   col_palette = NULL,
                                   text_color = "black") {

  # Extract graph and clusters from the output
  graph <- cluster_output$graph
  clusters <- cluster_output$clusters

  # Initialize a data frame for the table
  cluster_info <- list()

  # Helper function to get GO term names from GO IDs
  get_go_term_name <- function(go_id) {
    go_term_name <- tryCatch({
      AnnotationDbi::Term(go_id)
    }, error = function(e) {
      return(go_id)  # If the GO term is not found, return the GO ID itself
    })
    return(go_term_name)
  }

  if (method_type == "similarity") {
    if (is.null(similarity_matrix)) {
      stop("A GO semantic similarity object is required for similarity-based clustering.")
    }

    # For each cluster, determine the most representative pathway and list GO IDs
    for (cluster_id in unique(clusters)) {
      go_ids <- names(which(clusters == cluster_id))

      # Calculate the most representative pathway
      # Representative pathway: the term with the highest average similarity to other terms in the cluster
      sub_matrix <- similarity_matrix[go_ids, go_ids, drop = FALSE] # Subset similarity matrix
      avg_similarity <- rowMeans(sub_matrix)
      representative_pathway_id <- names(which.max(avg_similarity))

      # Get the full name of the representative pathway
      representative_pathway_name <- get_go_term_name(representative_pathway_id)

      # Check for NA values before proceeding
      if (is.na(representative_pathway_name) || length(go_ids) == 0) {
        next  # Skip this cluster if any values are NA
      }

      # Store information in the list
      cluster_info[[paste("Cluster", cluster_id)]] <- data.frame(
        Cluster = paste("Cluster", cluster_id),
        "Representative Pathway" = representative_pathway_name,
        "GO IDs" = paste(go_ids, collapse = ", "),
        stringsAsFactors = FALSE
      )
    }

    # Combine all cluster info into a single data frame
    cluster_df <- do.call(rbind, cluster_info)

  } else if (method_type == "network") {
    # For network-based clustering, simply list GO IDs for each cluster
    for (cluster_id in unique(clusters)) {
      go_ids <- names(which(clusters == cluster_id))

      # Check for NA values before proceeding
      if (length(go_ids) == 0) {
        next  # Skip this cluster if there are no GO IDs
      }

      # Store information in the list
      cluster_info[[paste("Cluster", cluster_id)]] <- data.frame(
        Cluster = paste("Cluster", cluster_id),
        "GO IDs" = paste(go_ids, collapse = ", "),
        stringsAsFactors = FALSE
      )
    }

    # Combine all cluster info into a single data frame
    cluster_df <- do.call(rbind, cluster_info)
  } else {
    stop("Invalid method_type. Choose either 'similarity' or 'network'.")
  }

  # Generate a table image
  # Use kableExtra for better table formatting
  if (!requireNamespace("kableExtra", quietly = TRUE)) {
    stop("The kableExtra package is required for table visualization. Install it using install.packages('kableExtra').")
  }
  library(kableExtra)

  # Add colors to rows based on clusters
  # Visualize clusters (assign colors)
  if(is.null(col_palette)){
    cluster_colors <- visualize_clusters(graph, col_palette)
    names(cluster_colors) <- paste("Cluster ", names(cluster_colors), sep = "")
  } else {
    cluster_colors <- col_palette
  }

  cluster_df %>%
    kbl(format = "html", escape = FALSE, row.names = FALSE) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = FALSE) %>%
    row_spec(0, bold = TRUE) %>%
    column_spec(1, color = text_color,
                background = cluster_colors[1:nrow(cluster_df)])
}
