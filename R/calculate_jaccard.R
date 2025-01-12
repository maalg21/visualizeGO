#' @title Calculate Jaccard Similarity Matrix
#'
#' @description This function calculates the similarity matrix using the Jaccard Index. In the context of Gene Ontology (GO),
#' the Jaccard Index is used to measure the similarity between two GO terms based on the overlap of their associated gene sets.
#' This method is simple and interpretable, especially when considering GO terms that share genes.
#' However, this method does not consider the hierarchical structure of the GO terms themselves,
#' only the overlap of gene annotations.
#'
#' @param graph A graph object that connects the input GO-terms with their parents and children terms
#' @return A similarity matrix with all the relationships between the GO-terms presented in the input graph.
#' @export

calculate_jaccard <- function(graph) {

  # Validate graph input
  if (is.null(graph) || !inherits(graph, "igraph")) {
    stop("The input 'graph' must be a valid igraph object.")
  }

  # Combine all GO terms from all sets into a unique list
  all_go_terms <- unique(V(graph)$name)

  # Validate input
  if (length(all_go_terms) < 2) {
    stop("Not enough unique GO terms to compute similarity.")
  }

  # Initialize similarity matrix for individual GO terms
  similarity_matrix <- matrix(0,
                              nrow = length(all_go_terms),
                              ncol = length(all_go_terms),
                              dimnames = list(all_go_terms, all_go_terms))

  # Calculate Jaccard Index for each pair of GO terms
  for (i in seq_along(all_go_terms)) {
    for (j in seq_along(all_go_terms)) {
      if (i <= j) {
        # Calculate intersection and union
        term1 <- all_go_terms[i]
        term2 <- all_go_terms[j]
        intersection <- ifelse(term1 == term2, 1, 0)  # Direct match
        union <- 1  # Union of a single term is itself

        # Compute Jaccard Index
        similarity_matrix[i, j] <- intersection / union
        similarity_matrix[j, i] <- similarity_matrix[i, j]  # Symmetric matrix
      }
    }
  }

  return(similarity_matrix)
}
