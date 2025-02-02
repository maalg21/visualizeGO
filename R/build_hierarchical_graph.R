#' @title Build Hierarchical Graph
#'
#' @description This function links all parents and children in the list of input GO terms and creates a graph, establishing hierarchical relationships between the GO terms.
#'
#' @param go_list1 List of GO-terms of interest.
#' @param go_list2 If enabled, second list of GO-terms to compare with the first one.
#' @param go_sim_object GO-terms similarity relationships save as environment.
#' @param nb_lists Allows to select whenever you are gonna use one ("single") or two ("double") lists.
#'
#' @return (1) A graph object that links the input GO-terms with their ancestors and childrens, (2).
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
      go_list1, go_list2,
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
      return(NULL)
    }
  }))

  # If no edges are found, include isolated terms in the graph
  isolated_terms <- setdiff(all_terms, unique(edges$to))
  if (length(isolated_terms) > 0) {
    isolated_edges <- data.frame(from = isolated_terms, to = isolated_terms)
    edges <- rbind(edges, isolated_edges)
  }

  # Remove self-loops
  edges <- edges[edges$from != edges$to, ]

  # Create graph object
  if (nrow(edges) > 0){
    graph <- graph_from_data_frame(edges, directed = TRUE)

    # Obtain input terms
    if(nb_lists == "single"){
      input_terms <- unique(go_list1)
    } else if (nb_lists == "double"){
      input_terms <- unique(c(go_list1, go_list2))
    }

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

  # Add the origin for each node
  if(nb_lists == "double"){
    V(graph)$list <- case_when(
      V(graph)$name %in% go_list1 ~ "List1",
      V(graph)$name %in% go_list2 ~ "List2",
      TRUE ~ NA)
  } else {
    V(graph)$list <- ifelse(V(graph)$name %in% go_list1, "List1", NA)
  }

  return(graph)
}
