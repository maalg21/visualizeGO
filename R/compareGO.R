#' @title Compare lists of GO-terms clusters
#'
#' @description This function aims to compare clusters of GO-terms that have been previously detected by semantic similarity. The comparison is performed by one of the similarity methods such as Resnik, Lin and Wang.
#'
#' @param comparison What kind of comparison is going to be made, whether "GO" or "cluster".
#' @param list1 A named vector that indicates the cluster to which each GO-term belongs,
#' i.e.: Variable clusters from the output of the cluster_go_terms function.
#' @param list2 Another set of GO-terms to compare with.
#' @param ontology Gene Ontology category to use (could be "BP" for Biological Process, "CC" for Cellular Component or "MF" for "Molecular Function").
#' @param OrgDb Organism to use as reference to obtain the GO-terms similarities. GOSemSimDATA object. Default = "org.Hs.eg.db" (human)
#' @param method One of "Resnik", "Lin" and "Wang" methods.
#' @param combine One of "max", "avg", "rcmax", "BMA" methods, for combining semantic similarity scores of multiple GO terms associated with protein or multiple proteins assiciated with protein cluster.
#' @param plot Option to display or not the graph. Default = TRUE
#' @param main Title of the HeatMap.
#' @param low Colour of the minimum value of the matrix
#' @param high Colour of the highest value of the matrix
#' @param labs Name of the X and Y axes.
#' @param values Option to display the similarity values for comparison face in the grid. Default = TRUE
#' @param cex Font size of the inside of the matrix in the plot.
#' @param cex_axis Font size of the axis title in the plot.
#'
#' @return Semantic similarity matrix of GO-terms lists or lists of GO-terms clusters.
#' If selected, the matrix is displayed in graph format.
#' @export

compareGO <- function(comparison = c("GO", "cluster"),
                      list1, list2,
                      ontology = c("BP", "CC", "MF"),
                      OrgDb = "org.Hs.eg.db",
                      method = c("Resnik", "Lin", "Wang"),
                      combine = c("avg", "BMA", "max", "rcmax"),
                      plot = T, low = "white", high = "red3",
                      labs = NULL, cex = 3, cex_axis = 10){

  # Load the OrgDb package dynamically ----
  if (!requireNamespace(OrgDb, quietly = TRUE)) {
    stop(paste("Package", OrgDb, "is required but not installed."))
  }
  library(OrgDb, character.only = TRUE)

  if(is.null(OrgDb)){
    stop("OrgDb object is required.")
  }
  cat("Getting GOSemSimData...\n")

  # Check validity of the method ----
  if (!method %in% c("Resnik", "Lin", "Wang")) {
    stop("Invalid method. Choose from 'Resnik', 'Lin', or 'Wang'.")
  }

  # Two lists for comparison ----
  if(is.null(list1) | is.null(list2)){
    stop("Two lists of GO-terms or of clustering results are needed for comparison based on semantic similarity.")
  }

  if (comparison == "cluster") {
    if (is.null(names(list1)) && all(names(list1) != "")) {
      stop("list1 is not a named vector.")
    }

    if (is.null(names(list2)) && all(names(list2) != "")) {
      stop("list2 is not a named vector.")
    }

    if (!combine %in% c("avg", "BMA", "max", "rcmax")){
      stop("One of 'max', 'avg', 'rcmax', 'BMA' methods")
    }
  }

  # Create a GO similarity object ----
  if(method != "Wang"){
    library(GOSemSim)
    go_data <- GOSemSim::godata(OrgDb = OrgDb, ont = ontology, computeIC = T)
  } else {
    go_data <- GOSemSim::godata(ont = ontology,  OrgDb = OrgDb, computeIC = FALSE)
  }

  # GO-terms comparison ----
  if(comparison == "GO"){
    # Initialize the similarity matrix
    semantic_similarity <- matrix(0, nrow = length(list1),
                                ncol = length(list2),
                                dimnames = list(list1, list2))

    # Calculate pairwise similarities
    for (i in seq_along(list1)) {
      for (j in seq_along(list2)) {
        semantic_similarity[i, j] <- GOSemSim::goSim(list1[i], list2[j],
                                                   semData = go_data,
                                                   measure = method)
      }
    }
  }

  # Clusters comparison ----
  if (comparison == "cluster") {
    m <- length(unique(list1))
    n <- length(unique(list2))

    semantic_similarity <- matrix(nrow = m, ncol = n)
    rownames(semantic_similarity) <- paste("Cluster ",
                                           1:length(unique(list1)),
                                           sep = "")
    colnames(semantic_similarity) <- paste("Cluster ",
                                           1:length(unique(list2)),
                                           sep = "")
    for(m in 1:m){
      for(n in 1:n){
        value <- GOSemSim::mgoSim(GO1 = names(list1[list1 == m]),
                                  GO2 = names(list2[list2 == n]),
                                  semData = go_data,
                                  measure = method,
                                  combine = combine)
        semantic_similarity[m, n] <- value
      }
    }
  }

  # Heatmap ----
  if(isTRUE(plot)){
    library(ggplot2)
    library(dplyr)
    library(tidyverse)
    p <- ggplot(data = semantic_similarity %>%
                  reshape2::melt(value.name = "SemanticSimilarity"),
                mapping = aes(x = Var1,
                              y = factor(Var2, levels = rev(unique(Var2))),
                              fill = SemanticSimilarity)) +
      geom_tile() + geom_text(aes(label = round(SemanticSimilarity, 2)),
                              size = cex, color = "black") +
      scale_fill_gradient(name = paste("GO: ", ontology,
                                       "\nSemantic Similarity\n(",
                                       method, ")", sep = ""),
                          low = low, high = high) +
      labs(x = labs[1], y = labs[2]) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = cex_axis),
            axis.text.y = element_text(color = "black", size = cex_axis),
            legend.title = element_text(face = "bold", hjust = .5),
            axis.title = element_text(face = "bold", size = 13),
            plot.title = element_text(face = "bold", size = 15, hjust = .5))
    return(p)
  } else {
    return(semantic_similarity)
  }
}
