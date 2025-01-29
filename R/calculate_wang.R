#' Calculate Wang Similarity Matrix
#'
#' This function calculates the similarity matrix using the Wang Method. The Wang Method is a
#' semantic similarity measure used to calculate the similarity between two Gene Ontology (GO)
#' terms based on the information content (IC) of their common ancestor and the terms themselves.
#' It incorporates both the GO hierarchy and the statistical properties of the terms to provide a
#' similarity score. The advantages of this method are that like other methods,
#' it combines ontology structure and statistical significance,
#' making it more accurate in capturing semantic similarity. However, requires accurate computation of
#' information content, which may be difficult if there is sparse data or limited term annotations.
#'
#' @param graph A graph object that connects the input GO-terms with their parents and children terms
#' @param ontology Gene Ontology category to use (could be "BP" for Biological Process, "CC" for Cellular Component or "MF" for "Molecular Function").
#' @param OrgDb Organism to use as reference to obtain the GO-terms similarities. Default = "org.Hs.eg.db"
#' @return A similarity matrix with all the relationships between the GO-terms presented in the input graph.
#' @references Wang et al., (2007) A new method to measure the semantic similarity of GO terms. *Bioinformatics*, 23:10, 1274–1281, \href{https://doi.org/10.1093/bioinformatics/btm087}{10.1093/bioinformatics/btm087}
#' @export

calculate_wang <- function(graph, ontology = c("BP", "CC", "MF"),
                                      OrgDb = "org.Hs.eg.db") {
  # Load required libraries
  if (!requireNamespace("GOSemSim", quietly = TRUE)) {
    stop("Please install the 'GOSemSim' package.")
  }

  # Load GO Semantic Similarity data
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

  sem_data <- godata(ont = ontology, OrgDb = OrgDb, computeIC = FALSE)

  # Initialize similarity matrix for all GO terms
  similarity_matrix <- matrix(0,
                              nrow = length(all_go_terms),
                              ncol = length(all_go_terms),
                              dimnames = list(all_go_terms, all_go_terms))

  # Calculate Wang Similarity for each pair of sets
  for (i in seq_along(all_go_terms)) {
    for (j in seq_along(all_go_terms)) {
      if (i <= j) {
        # Compute Wang Similarity between two GO-term sets
        sim <- mgoSim(
          all_go_terms[i], all_go_terms[j],
          semData = sem_data, measure = "Wang", combine = "avg")

        similarity_matrix[i, j] <- sim
        similarity_matrix[j, i] <- sim  # Symmetric matrix
      }
    }
  }

  return(similarity_matrix)
}
