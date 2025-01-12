#' Retain Ancestors
#'
#' This function reload the parent input GO-terms and select those direct ancestors that in the filtering before could be removed.
#'
#' @param graph A graph object from the "build_hierarchical_graph" or "expand_graph" functions
#' @param go_list1 List of GO-terms of interest.
#' @param go_list2 If enabled, second list of GO-terms to compare with the first one.
#' @param nb_lists Allows to select whenever you are gonna use one ("single") or two ("double") lists
#' @return A graph object that updated with all the direct ancestors of the input GO-terms.
#' @export

retain_ancestors_above_input_terms <- function(graph, go_list1, go_list2 = NULL,
                                               nb_lists = c("single", "double")) {

  # Match argument for nb_lists
  nb_lists <- match.arg(nb_lists)

  # Handle single or double list cases
  if (nb_lists == "single") {
    relevant_go_terms <- unique(go_list1)
  } else if (nb_lists == "double") {
    if (is.null(go_list2)) {
      stop("When nb_lists is 'double', both go_list1 and go_list2 must be provided.")
    }
    relevant_go_terms <- unique(c(go_list1, go_list2))
  } else {
    stop("Invalid value for nb_lists. Use 'single' or 'double'.")
  }

  # Identify the ancestors (parent GO-terms) of the input GO-terms
  ancestors <- unique(unlist(lapply(relevant_go_terms, function(term) {
    # Extracting neighbors with 'mode = "in"' for ancestors (parent terms)
    ancestor_indices <- neighbors(graph, term, mode = "in")

    # Convert numeric indices to GO-term IDs (names)
    ancestor_terms <- V(graph)$name[ancestor_indices]

    return(ancestor_terms)  # Return the GO-term IDs of the ancestors
  })))

  # Include the input terms themselves and their ancestors
  terms_to_keep <- unique(c(relevant_go_terms, ancestors))

  # Now, to ensure the graph only contains these terms and their relationships,
  # we will extract all the nodes that are related (either as input terms or ancestors)
  filtered_graph <- induced_subgraph(graph, vids = V(graph)$name %in% terms_to_keep)

  return(filtered_graph)
}
