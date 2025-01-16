#' @title Compare lists of GO-terms clusters
#'
#' @description asdfghjlñ
#'
#' @param go_list1 List of GO-terms of interest.
#' @param go_list2 If enabled, second list of GO-terms to compare with the first one.
#' @param go_sim_object GO-terms similarity relationships save as environment.
#' @param nb_lists Allows to select whenever you are gonna use one ("single") or two ("double") lists
#' @return A graph object that links the input GO-terms with their ancestors and childrens.
#' @export

compare_clusters <- function(cluster_list1, cluster_list2,
                             ontology = c("BP", "CC", "MF"),
                             OrgDb = "org.Hs.eg.db"){
  hsGO <- godata(OrgDb = 'org.Hs.eg.db', ont = "BP")  # Biological Process ontology

  # List of GO-terms clusters
  Perirenal <- split(names(clusterPF$clusters),
                     paste("Cluster ", clusterPF$clusters, sep = ""))
  Tail <- split(names(clusterTF$clusters),
                paste("Cluster ", clusterTF$clusters, sep = ""))

  Pe <- length(Perirenal)
  Ta <- length(Tail)

  semantic_similarity <- matrix(nrow = Pe, ncol = Ta)
  rownames(semantic_similarity) <- paste("Cluster ", 1:Pe, sep = "")
  colnames(semantic_similarity) <- paste("Cluster ", 1:Ta, sep = "")
  for(Pe in 1:Pe){
    for(Ta in 1:Ta){
      value <- mgoSim(GO1 = unlist(Perirenal[Pe]),
                      GO2 = unlist(Tail[Ta]),
                      semData = hsGO,
                      measure = "Wang",
                      combine = "BMA")
      semantic_similarity[Pe, Ta] <- value
    }
  }

  library(ggplot2)
  library(dplyr)
  library(tidyverse)
  ggplot(data = semantic_similarity %>%
           reshape2::melt(value.name = "SemanticSimilarity") %>%
           dplyr::rename("Perirenal" = Var1,
                         "Tail" = Var2),
         mapping = aes(x = Perirenal,
                       y = factor(Tail, levels = rev(unique(Tail))),
                       fill = SemanticSimilarity)) +
    geom_tile() + geom_text(aes(label = SemanticSimilarity), size = 3) +
    scale_fill_gradient(name = "Semantic\nSimilarity",
                        low = "white", high = "red3") +
    labs(x = "Perirenal Clusters", y = "Tail Clusters") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"))

}
