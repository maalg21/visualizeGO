#' @title Compare lists of GO-terms clusters
#'
#' @description This function aims to compare clusters of GO-terms that have been previously detected by semantic similarity. The comparison is performed by one of the similarity methods such as Resnik, Lin and Wang.
#'
#' @param cluster_list1 A set of GO-terms.
#' @param cluster_list2 Another set of GO-terms.
#' @param ontology Gene Ontology category to use (could be "BP" for Biological Process, "CC" for Cellular Component or "MF" for "Molecular Function").
#' @param OrgDb Organism to use as reference to obtain the GO-terms similarities. GOSemSimDATA object. Default = "org.Hs.eg.db" (human)
#' @param method One of "Resnik", "Lin" and "Wang" methods.
#' @param combine One of "max", "avg", "rcmax", "BMA" methods, for combining semantic similarity scores of multiple GO terms associated with protein or multiple proteins assiciated with protein cluster.
#' @param plot Option to display or not the graph. Default = TRUE
#' @param low Colour of the minimum value of the matrix
#' @param high Colour of the highest value of the matrix
#' @param labs Name of the X and Y axes.
#' @param cex Font size of the inside of the matrix in the plot.
#'
#' @return Semantic similarity matrix of GO-terms clusters. If selected, the matrix is displayed in graph format.
#' @export

compare_clusters <- function(cluster_list1, cluster_list2,
                             ontology = c("BP", "CC", "MF"),
                             OrgDb = "org.Hs.eg.db",
                             method = c("Resnik", "Lin", "Wang"),
                             combine = c("avg", "BMA", "max", "rcmax"),
                             plot = T, low = "white", high = "red3",
                             labs = NULL, cex = 3){

  # Step 1. Obtaining the GOSemSimDATA ----
  if(is.null(OrgDb)){
    stop("OrgDb object is required.")
  }
  cat("Getting GOSemSimData...\n")
  GOData <- godata(annoDb = OrgDb, ont = ontology)

  # Step 2. List of GO-terms clusters ----
  if(is.null(cluster_list1) | is.null(cluster_list2)){
    stop("Two clustering results are needed through semantic similarity")
  }

  if (!"clusters" %in% names(cluster_list1) || !"clusters" %in% names(cluster_list2)) {
    stop("The clusters input must contain 'clusters' components.")
  }

  List1 <- split(names(cluster_list1$clusters),
                     paste("Cluster ", cluster_list1$clusters, sep = ""))
  List2 <- split(names(cluster_list2$clusters),
                paste("Cluster ", cluster_list2$clusters, sep = ""))

  m <- length(List1)
  n <- length(List2)

  semantic_similarity <- matrix(nrow = m, ncol = n)
  rownames(semantic_similarity) <- paste("Cluster ", 1:m, sep = "")
  colnames(semantic_similarity) <- paste("Cluster ", 1:n, sep = "")
  for(m in 1:m){
    for(n in 1:n){
      value <- mgoSim(GO1 = unlist(Perirenal[m]),
                      GO2 = unlist(Tail[n]),
                      semData = GOData,
                      measure = method,
                      combine = combine)
      semantic_similarity[m, n] <- value
    }
  }

  if(isTRUE(plot)){
    library(ggplot2)
    library(dplyr)
    library(tidyverse)
    p <- ggplot(data = semantic_similarity %>%
             reshape2::melt(value.name = "SemanticSimilarity"),
           mapping = aes(x = Var1,
                         y = factor(Var2, levels = rev(unique(Var2))),
                         fill = SemanticSimilarity)) +
      geom_tile() + geom_text(aes(label = SemanticSimilarity), size = cex) +
      scale_fill_gradient(name = "Semantic\nSimilarity",
                          low = low, high = high) +
      labs(x = labs[1], y = labs[2]) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.title = element_text(face = "bold"),
            axis.title = element_text(face = "bold"))
    return(p)
  } else {
    return(semantic_similarity)
  }
}
