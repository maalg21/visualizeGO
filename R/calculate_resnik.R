#' Calculate Resnik Similarity Matrix
#'
#' This function calculates the similarity matrix using the Resnik Method.
#' The Resnik Method is a technique used to calculate the semantic similarity between two Gene Ontology (GO)
#' terms based on their shared information content (IC). The core idea is that the more specific (or lower)
#' the most specific common ancestor of the two GO terms, the higher their semantic similarity.
#' It's simple and intuitive, based on the notion that common ancestors reflect shared biological meanings.
#' However, relies heavily on the structure of the GO ontology and can be influenced by how well the dataset represents the real biological phenomena.
#'
#' @param graph A graph object that connects the input GO-terms with their parents and children terms
#' @param ontology Gene Ontology category to use (could be "BP" for Biological Process, "CC" for Cellular Component or "MF" for "Molecular Function").
#' @param OrgDb Organism to use as reference to obtain the GO-terms similarities. Default = "org.Hs.eg.db"
#' @return A similarity matrix with all the relationships between the GO-terms presented in the input graph.
#' @references Resnik (1995) Using Information Content to Evaluate Semantic Similarity in a Taxonomy. *Proceedings of the 14th International Joint Conference on Artificial Intelligence*, \href{https://doi.org/10.48550/arXiv.cmp-lg/9511007}{10.48550/arXiv.cmp-lg/9511007}
#' @export

calculate_resnik <- function(graph, ontology = c("BP", "CC", "MF"),
                             OrgDb = "org.Hs.eg.db") {
  # Load necessary library
  if (!requireNamespace("GOSemSim", quietly = TRUE)) {
    stop("The GOSemSim package is required. Install it using install.packages('GOSemSim').")
  }

  library(GOSemSim)

  # Load the appropriate GO data
  if (!requireNamespace(OrgDb, quietly = TRUE)) {
    stop(paste("Please install the", OrgDb, "package to proceed."))
  }

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

  # Prepare the GO data
  sem_data <- godata(OrgDb = OrgDb, ont = ontology, computeIC = TRUE)

  # Initialize similarity matrix for individual GO terms
  similarity_matrix <- matrix(0,
                              nrow = length(all_go_terms),
                              ncol = length(all_go_terms),
                              dimnames = list(all_go_terms, all_go_terms))

  # Compute pairwise Resnik similarity
  for (i in seq_along(all_go_terms)) {
    for (j in seq_along(all_go_terms)) {
      if (i <= j) {
        # Calculate Resnik similarity between two individual GO terms
        sim <- goSim(all_go_terms[i], all_go_terms[j], semData = sem_data, measure = "Resnik")
        similarity_matrix[i, j] <- sim
        similarity_matrix[j, i] <- sim  # Symmetric matrix
      }
    }
  }

  return(similarity_matrix)
}
