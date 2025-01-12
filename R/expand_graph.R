#' Expand Graph
#'
#' This function filters the input graph to remove all disconnected GO-terms
#'
#' @param graph A graph object from the "build_hierarchical_graph" function
#' @param go_list1 List of GO-terms of interest.
#' @param go_list2 If enabled, second list of GO-terms to compare with the first one.
#' @param nb_lists Allows to select whenever you are gonna use one ("single") or two ("double") lists
#' @return A graph object that only manteins the GO-terms connected with the input GO-terms.
#' @export

expand_graph <- function(graph, go_list1, go_list2 = NULL,
                         nb_lists = c("single", "double")) {

  # Match argument for nb_lists
  nb_lists <- match.arg(nb_lists)

  # Handle single or double list cases
  if (nb_lists == "single") {
    input_terms <- unique(go_list1)
  } else if (nb_lists == "double") {
    if (is.null(go_list2)) {
      stop("When nb_lists is 'double', both go_list1 and go_list2 must be provided.")
    }
    input_terms <- unique(c(go_list1, go_list2))
  } else {
    stop("Invalid value for nb_lists. Use 'single' or 'double'.")
  }

  # Check if all input terms are in the graph
  if (!all(input_terms %in% V(graph)$name)) {
    missing_terms <- input_terms[!(input_terms %in% V(graph)$name)]
    warning("The following terms are missing from the graph: ", paste(missing_terms, collapse = ", "))
    # Optionally add missing terms back to the graph
    graph <- add_vertices(graph, length(missing_terms), name = missing_terms)
  }

  # Get neighbors for the connected terms (both ancestors and descendants)
  all_ancestors <- unique(unlist(lapply(input_terms, function(term) {
    # Extracting neighbors with 'mode = "in"' for ancestors
    ancestors <- neighbors(graph, term, mode = "in")
    valid_ancestors <- V(graph)$name[ancestors]
    return(valid_ancestors)  # Get names of the ancestor nodes
  })))

  all_descendants <- unique(unlist(lapply(input_terms, function(term) {
    # Extracting neighbors with 'mode = "out"' for descendants
    descendants <- neighbors(graph, term, mode = "out")
    valid_descendants <- V(graph)$name[descendants]
    return(valid_descendants)  # Get names of the descendant nodes
  })))

  # Combine the relevant terms
  all_relevant_terms <- unique(c(input_terms, all_ancestors, all_descendants))

  # Filter the graph to include only the relevant terms and their relationships
  filtered_graph <- induced_subgraph(graph,
                                     vids = V(graph)[V(graph)$name %in%
                                                       all_relevant_terms])

  return(filtered_graph)
}
