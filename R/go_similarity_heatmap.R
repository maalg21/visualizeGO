#' GO Similarity HeatMap
#'
#' This function compares the similarity between two lists of GO-terms through any of the chosen methods.
#'
#' @param go_list1 First list of GO-terms to compare.
#' @param go_list2 Second list of GO-terms to compare.
#' @param ont Gene Ontology category to use (could be "BP" for Biological Process, "CC" for Cellular Component or "MF" for "Molecular Function").
#' @param orgdb Organism to use as reference to obtain the GO-terms similarities. Default = "org.Hs.eg.db"
#' @param xlab Name of the X-axis.
#' @param ylab Name of the Y-axis.
#' @param cex Size of the labels.
#' @param main Title of the HeatMap.
#' @return A HeatMap is obtained whose intensity is given by the degree of similarity between the GO-terms being compared.
#' @export

go_similarity_heatmap <- function(go_list1, go_list2,
                                  ont = c("BP", "MF", "CC"),
                                  method = c("Jaccard", "Resnik", "Lin", "Wang"),
                                  orgdb = "org.Hs.eg.db",
                                  xlab = NULL, ylab = NULL,
                                  cex = 10, main = NULL) {

  # Check validity of the method
  if (!method %in% c("Jaccard", "Resnik", "Lin", "Wang")) {
    stop("Invalid method. Choose from 'Jaccard', 'Resnik', 'Lin', or 'Wang'.")
  }

  # Load the OrgDb package dynamically
  if (!requireNamespace(orgdb, quietly = TRUE)) {
    stop(paste("Package", orgdb, "is required but not installed."))
  }
  library(orgdb, character.only = TRUE)

  if(method != "Jaccard"){
    if(method != "Wang"){
      # Create a GO similarity object
      go_data <- GOSemSim::godata(orgdb, ont = ont,
                                  computeIC = method %in% c("Resnik", "Lin"))
    } else {
      go_data <- GOSemSim::godata(ont = ontology, OrgDb = orgdb, computeIC = FALSE)
    }

    # Initialize the similarity matrix
    similarity_matrix <- matrix(0, nrow = length(go_list1),
                                ncol = length(go_list2),
                                dimnames = list(go_list1, go_list2))

    # Calculate pairwise similarities
    for (i in seq_along(go_list1)) {
      for (j in seq_along(go_list2)) {
        similarity_matrix[i, j] <- goSim(go_list1[i], go_list2[j],
                                         semData = go_data,
                                         measure = method)
      }
  }
  } else {
    for (i in seq_along(go_list1)) {
      for (j in seq_along(go_list2)) {
        if (i <= j) {
          # Calculate intersection and union
          term1 <- go_list1[i]
          term2 <- go_list2[j]
          intersection <- ifelse(term1 == term2, 1, 0)  # Direct match
          union <- 1  # Union of a single term is itself

          # Compute Jaccard Index
          similarity_matrix[i, j] <- intersection / union
          similarity_matrix[j, i] <- similarity_matrix[i, j]  # Symmetric matrix
        }
      }
    }
  }

  # Generate the heatmap
  heatmap_plot <- ggplot(data = as.data.frame(similarity_matrix) %>%
                           tibble::rownames_to_column(var = "GO_ID1") %>%
                           reshape2::melt(value.name = "SemanticSimilarity") %>%
                           dplyr::rename("GO_ID2" = variable)) +
    geom_tile(mapping = aes(x = GO_ID1, y = GO_ID2, fill = SemanticSimilarity)) +
    scale_fill_gradient(name = paste("GO: ", ont,
                                     "\nSemantic Similarity\n(",
                                     method, ")", sep = ""),
                        low = "white", high = "red") +
    labs(x = xlab, y = ylab) + theme_minimal() +
    ggtitle(main) +
    theme(axis.text.x = element_text(angle = 45, color = "black",
                                     hjust = 1, size = cex),
          axis.text.y = element_text(color = "black", size = cex),
          legend.title = element_text(hjust = .5, face = "bold"),
          axis.title = element_text(face = "bold", size = 13),
          plot.title = element_text(face = "bold", size = 15, hjust = .5))

  print(heatmap_plot)

  return(similarity_matrix)
}
