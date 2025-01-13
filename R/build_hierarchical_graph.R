#' @title Build Hierarchical Graph
#'
#' @description This function link all the parent and child of the input GO-terms list and make a graph
#'
#' @param go_list1 List of GO-terms of interest.
#' @param go_list2 If enabled, second list of GO-terms to compare with the first one.
#' @param go_sim_object GO-terms similarity relationships save as environment.
#' @param nb_lists Allows to select whenever you are gonna use one ("single") or two ("double") lists
#' @return A graph object that links the input GO-terms with their ancestors and childrens.
#' @export

build_hierarchical_graph <- function(go_list1, go_list2 = NULL, go_sim_object = NULL,
                                     nb_lists = c("single", "double")) {
  # Match argument for nb_lists
  nb_lists <- match.arg(nb_lists)

  # Ensure igraph is loaded
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("The 'igraph' package is required but not installed. Please install it using install.packages('igraph').")
  }
  library(igraph)

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

  # Handle single or double list cases
  if (nb_lists == "single") {
    all_terms <- unique(c(
      go_list1,
      get_all_ancestors(go_list1, go_sim_object)
    ))
  } else if (nb_lists == "double") {
    if (is.null(go_list2)) {
      stop("When nb_lists is 'double', both go_list1 and go_list2 must be provided.")
    }
    all_terms <- unique(c(
      go_list1,
      go_list2,
      get_all_ancestors(go_list1, go_sim_object),
      get_all_ancestors(go_list2, go_sim_object)
    ))
  } else {
    stop("Invalid value for nb_lists. Use 'single' or 'double'.")
  }

  # Build edges between terms and their children
  edges <- do.call(rbind, lapply(all_terms, function(term) {
    children <- get_children(term, go_sim_object)
    if (!is.null(children)) {
      return(data.frame(from = term, to = children))
    } else {
      return(NULL)  # Skip if no children exist
    }
  }))

  # If no edges are found, include isolated terms in the graph
  isolated_terms <- setdiff(all_terms, unique(edges$to))
  if (length(isolated_terms) > 0) {
    isolated_edges <- data.frame(from = isolated_terms, to = isolated_terms)
    edges <- rbind(edges, isolated_edges)
  }

  # Remove self-loops (edges where from == to)
  edges <- edges[edges$from != edges$to, ]

  # Create graph object if edges exist
  if (nrow(edges) > 0) {
    graph <- graph_from_data_frame(edges, directed = TRUE)
    return(graph)
  } else {
    stop("No edges found for the provided GO terms.")
  }
}
