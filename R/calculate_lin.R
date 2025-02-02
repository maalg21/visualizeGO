#' Calculate Lin Similarity Matrix
#'
#' This function calculates the similarity matrix using the Lin Method.
#' The Lin Method is a semantic similarity measure for calculating the similarity between two
#' Gene Ontology (GO) terms based on their Most Specific Common Ancestor (MSCA) and information content (IC).
#' It incorporates both the structure of the GO hierarchy and the frequency of terms in the ontology
#' to quantify similarity. This method accounts for both the structure of the ontology (using the MSCA)
#' and the frequency of terms, making it a balanced method for calculating similarity. However,
#' relies on the availability of accurate information content values and can be sensitive to the way the ontology is structured and annotated.
#'
#' @param graph A graph object that connects the input GO-terms with their parents and children terms
#' @param ontology Gene Ontology category to use (could be "BP" for Biological Process, "CC" for Cellular Component or "MF" for "Molecular Function").
#' @param OrgDb Organism to use as reference to obtain the GO-terms similarities. Default = "org.Hs.eg.db"
#' @return A similarity matrix with all the relationships between the GO-terms presented in the input graph.
#' @references Lin (1998) An Information-Theoretic Definition of Similarity. *Proceedings of the Fifteenth International Conference on Machine Learning*
#' @export

calculate_lin <- function(graph, ontology = c("BP", "CC", "MF"),
                                     OrgDb = "org.Hs.eg.db") {
  # Load necessary library
  if (!requireNamespace("GOSemSim", quietly = TRUE)) {
    stop("The GOSemSim package is required. Install it using install.packages('GOSemSim').")
  }

  library(GOSemSim)

  # Validate inputs
  if (!is.vector(go_terms) || length(go_terms) < 2) {
    stop("Input must be a vector of at least two GO terms.")
  }

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

  # Create a semantic similarity dataset
  sem_data <- godata(OrgDb, ont = ontology, computeIC = TRUE)

  # Initialize similarity matrix for individual GO terms
  similarity_matrix <- matrix(0,
                              nrow = length(all_go_terms),
                              ncol = length(all_go_terms),
                              dimnames = list(all_go_terms, all_go_terms))

  # Compute pairwise Lin similarity
  for (i in seq_along(all_go_terms)) {
    for (j in seq_along(all_go_terms)) {
      if (i <= j) {
        # Calculate Lin similarity between two individual GO terms
        sim <- goSim(all_go_terms[i], all_go_terms[j], semData = sem_data, measure = "Lin")
        similarity_matrix[i, j] <- sim
        similarity_matrix[j, i] <- sim  # Symmetric matrix
      }
    }
  }

  return(similarity_matrix)
}
