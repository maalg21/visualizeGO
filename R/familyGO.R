#' @title Build Hierarchical Graph
#'
#' @description This function links all parents and children in the list of input GO terms and creates a graph, establishing hierarchical relationships between the GO terms.
#'
#' @param cluster Cluster object, obtained from clusterGO.
#' @param go_sim_object GO-terms similarity relationships save as environment.
#'
#' @return A igraph object that links the input GO-terms with their ancestors and childrens.
#' @export

familyGO <- function(cluster, go_sim_object = NULL){

  # Ensure igraph is loaded
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("The 'igraph' package is required but not installed. Please install it using install.packages('igraph').")
  }
  library(igraph)

  if (is.null(go_sim_object)) {
    file_name <- "go_term_database_all_ontologies.RData"
    file <- system.file("extdata", file_name, package = "visualizeGO")

    if (!file.exists(file)) {
      stop(paste("Error: The file", file_name, "does not exist in the directory or it's corrupted."))
    } else {
      message("File found: ", file_name)
      load(file)
      go_sim_object <- GOSimEnv
    }
  }

  if (is.null(cluster) ||
      !"clusters" %in% names(cluster) ||
      !"descriptions" %in% names(cluster) ||
      !"representative_pathways" %in% names(cluster)) {
    stop("Invalid cluster result provided.")
  }

  input_terms <- cluster$descriptions$GO_ID

  # Build edges between terms and their children
  edges <- do.call(rbind, lapply(input_terms, function(term) {
    children <- get_children(term, go_sim_object)
    if (!is.null(children)) {
      return(data.frame(from = term, to = children))
    } else {
      return(NULL)
    }
  }))

  # If no edges are found, include isolated terms in the graph
  isolated_terms <- setdiff(input_terms, unique(edges$to))
  if (length(isolated_terms) > 0) {
    isolated_edges <- data.frame(from = isolated_terms, to = isolated_terms)
    edges <- rbind(edges, isolated_edges)
  }

  # Remove self-loops
  edges <- edges[edges$from != edges$to, ]

  # Create graph object
  if (nrow(edges) > 0){
    graph <- graph_from_data_frame(edges, directed = TRUE)

    # Count parent occurrences for input terms
    parent_counts <- table(edges$from[edges$to %in% input_terms])

    # Count child occurrences for input terms
    child_counts <- table(edges$to[edges$from %in% input_terms])

    # Identify intermediate terms that are both parent and child of input terms
    intermediate_terms <- intersect(names(parent_counts),
                                    names(child_counts))

    important_parents <- names(parent_counts[parent_counts > 1])
    important_children <- names(child_counts[child_counts > 1])

    # Keep only relevant nodes
    relevant_terms <- unique(c(input_terms,
                               important_parents,
                               important_children,
                               intermediate_terms))

    # Ensure relevant terms exist in the graph before filtering
    relevant_terms <- relevant_terms[relevant_terms %in% V(graph)$name]

    if (length(relevant_terms) > 0) {
      graph <- induced_subgraph(graph, vids = V(graph)$name %in% relevant_terms)
    } else {
      warning("Filtering resulted in an empty graph. Returning the full graph instead.")
    }

    V(graph)$origin <- ifelse(V(graph)$name %in% input_terms, "input", "external")

  } else {
    stop("No edges found for the provided GO terms.")
  }

  # Add cluster information
  cluster_assignments <- as.data.frame(cluster$clusters) %>%
    tibble::rownames_to_column(var = "GO_ID") %>%
    dplyr::rename("Cluster" = `cluster$clusters`) %>%
    arrange(Cluster)

  cluster_assignments[cluster_assignments$GO_ID == g,]$Cluster

  V(graph)$Cluster <- ifelse(
    V(graph)$origin == "input",
    cluster_assignments$Cluster[match(V(graph)$name,
                                      cluster_assignments$GO_ID)],
    NA
  )

  return(graph)

}
