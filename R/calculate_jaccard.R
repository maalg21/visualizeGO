#' @title Calculate Jaccard Similarity Matrix
#'
#' @description This function calculates the similarity matrix using the Jaccard Index. In the context of Gene Ontology (GO),
#' the Jaccard Index is used to measure the similarity between two GO terms based on the overlap of their associated gene sets.
#' This method is simple and interpretable, especially when considering GO terms that share genes.
#' However, this method does not consider the hierarchical structure of the GO terms themselves,
#' only the overlap of gene annotations.
#'
#' @param input_terms A vector containing the input GO-terms IDs
#' @return A similarity matrix with all the relationships between the GO-terms calculated based on the Jaccard Index.
#' @export

calculate_jaccard <- function(input_terms) {

  # Validate input
  if (length(input_terms) < 2) {
    stop("Not enough unique GO terms to compute similarity.")
  }

  # Initialize similarity matrix for individual GO terms
  similarity_matrix <- matrix(0,
                              nrow = length(input_terms),
                              ncol = length(input_terms),
                              dimnames = list(input_terms, input_terms))

  # Calculate Jaccard Index for each pair of GO terms
  for (i in seq_along(input_terms)) {
    for (j in seq_along(input_terms)) {
      if (i <= j) {
        # Calculate intersection and union
        term1 <- input_terms[i]
        term2 <- input_terms[j]
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
