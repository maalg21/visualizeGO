#' Determine the Number of Clusters
#'
#' This function determines the number of clusters into which the GO-terms in the graph are divided based on the Elbow and Silhouette methods.
#'
#' @param similarity_matrix Only if the method_type selected was "similarity". Similarity matrix previously calculated that relates the GO-terms depending on their semantic similarity.
#' @param graph igraph object linking the input GO-terms.
#' @return A plot showing the optimal number of clusters based on the different methods.
#' @export

determine_nbclusters <- function(graph, similarity_matrix){

  if (!requireNamespace("factoextra", quietly = TRUE)) {
    stop("The factoextra package is required. Install it using install.packages('factoextra').")
  } else {
    library(factoextra)
  }

  distance_matrix <- 1 - similarity_matrix
  nb_terms <- length(V(graph)$name)

  # Elbow's method ----
  elbow <- fviz_nbclust(similarity_matrix, kmeans,
                        method = "wss",
                        diss = distance_matrix,
                        k.max = nb_terms-1)

  elbow_values <- elbow$data$y
  # Compute the first derivative (rate of change)
  diff <- diff(elbow_values)
  # Compute the second derivative (acceleration)
  diff2 <- diff(diff)
  # Find the knee/elbow point (where the second derivative is maximized)
  optimal_k_elbow <- which.max(abs(diff2)) + 1
  cat("Based on the Elbow method, the optimal number of clusters (k) is:", optimal_k_elbow, "\n")

  p1 <- elbow + theme_minimal() + ggtitle("The Elbow Method") +
    geom_line(color = "red", linewidth = 0.5, group = 1) +
    geom_point(color = "red") +
    geom_vline(xintercept = optimal_k_elbow, color = "red") +
    scale_x_discrete(breaks = as.character(seq(0, nb_terms, by = 20))) +
    labs(x = "Number of Clusters", y = "Total WSS", tag = "A") +
    theme(plot.title = element_text(face = "bold", hjust = .5, size = 13),
          axis.title = element_text(face = "bold"),
          plot.tag = element_text(face = "bold"))


  # Silhouette method ----
  sil <- fviz_nbclust(similarity_matrix, kmeans,
                      method = "silhouette",
                      diss = distance_matrix,
                      k.max = nb_terms-1)

  sil_scores <- sil$data$y
  optimal_k_sil <- which.max(sil_scores)
  cat("Based on the Silhouette method, the optimal number of clusters (k) is:", optimal_k_sil, "\n")

  p2 <- sil + theme_minimal() + ggtitle("The Silhouette Method") +
    geom_line(color = "red", linewidth = 0.5, group = 1) +
    geom_point(color = "red") +
    geom_vline(xintercept = optimal_k_sil, color = "red") +
    scale_x_discrete(breaks = as.character(seq(0, nb_terms, by = 20))) +
    labs(x = "Number of Clusters", y = "Average Silhouette Width", tag = "B") +
    theme(plot.title = element_text(face = "bold", hjust = .5, size = 13),
          axis.title = element_text(face = "bold"),
          plot.tag = element_text(face = "bold"))

  library(grid)
  library(gridExtra)
  grid.newpage()
  plot <- grid.arrange(p1, p2, ncol = 2)

  return(list(optimal_Silhouette = optimal_k_sil,
              optimal_Elbow = optimal_k_elbow, plot))
  print(plot)
}
