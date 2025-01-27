#' Scatter Plot of GO-terms
#'
#' Function to make a plot of GO terms as scattered points.
#' Distances between points represent the similarity between terms, and axes
#' are the first 2 components of applying a PCA to the similarity matrix.
#' Size of the point represents the number of genes the GO term contains.
#' Each color represents a cluster.
#'
#' @param similarity_matrix Squared matrix of the semantic similarity for each pair of terms calculated through any of the methods.
#' @param cluster Grouping of GO-terms in the different clusters.
#' @param title Title of the plot
#' @param colors Vector of colors for each of the represented clusters.
#' @param labels Enable to write the Clusters number in the plot
#'
#' @return A description of what the function returns (e.g., a numeric vector, data frame, list).
#'
#' @export
scatterGO <- function(similarity_matrix, cluster,
                      title = NULL, colors = NULL,
                      labels = T) {

  if (is.null(similarity_matrix)) {
    stop("A GO semantic similarity matrix is required for performing the plot.")
  }

  if (!"clusters" %in% names(clusters) || !"graph" %in% names(clusters)) {
    stop("The 'clusters' input must contain both 'graph' and 'clusters' components.")
  }

  distance_matrix <- 1 - similarity_matrix
  mds_result <- cmdscale(as.dist(distance_matrix), k = 2)
  pca_result <- prcomp(as.dist(distance_matrix))
  pca_scores <- pca_result$x
  pca_scores <- merge(pca_scores, by.x = "row.names",
                      data.frame("IDs" = names(cluster$clusters),
                                 "Cluster" = cluster$clusters), by.y = "IDs",
                      all = T) %>%
    tibble::column_to_rownames(var = "Row.names") %>%
    dplyr::mutate("Cluster" = paste("Cluster ", .$Cluster, sep = ""))

  library(ggplot2)
  if(is.null(colors)){
    colors <- generate_pastel_colors(n = cluster$nb_clusters)
  }

  if(isTRUE(labels)){
    library(dplyr)
    centroids <- pca_scores %>%
      group_by(Cluster) %>%
      dplyr::summarise(PC1 = mean(PC1), PC2 = mean(PC2))

    ggplot(pca_scores, aes(x = PC1, y = PC2,
                           color = Cluster)) +
      geom_point(size = 4, alpha = 0.8) +
      geom_label(data = centroids, aes(label = Cluster),
                 vjust = -0.5, hjust = 0.5) +
      scale_color_manual(values = colors) +
      geom_hline(yintercept = 0, colour = "black") +
      geom_vline(xintercept = 0, colour = "black") +
      theme_minimal() +
      theme(legend.position = "none",
            axis.title = element_text(size = 15, face = "bold", hjust = .5),
            plot.title = element_text(size = 20, face = "bold", hjust = .5)) +
      labs(title = title,
           x = "PC1",y = "PC2")
  } else {
    ggplot(pca_scores, aes(x = PC1, y = PC2,
                           color = Cluster)) +
      geom_point(size = 4, alpha = 0.8) +
      scale_color_manual(values = colors) +
      geom_hline(yintercept = 0, colour = "black") +
      geom_vline(xintercept = 0, colour = "black") +
      theme_minimal() +
      theme(legend.position = "none",
            axis.title = element_text(size = 15, face = "bold", hjust = .5),
            plot.title = element_text(size = 20, face = "bold", hjust = .5)) +
      labs(title = title,
           x = "PC1",y = "PC2")
  }
}
