#' Cluster GO Terms
#'
#' This function performs clustering of Gene Ontology (GO) terms based on either semantic similarity or network-based methods.
#'
#' @param method_type Choose between "similarity" or "network" to select the method to assign the clusters to each input GO term.
#' @param method Once selected the type of assignment, it's the method within each type selected to obtain the clusters.
#' @param similarity_matrix Only if the method_type selected was "similarity". Similarity matrix previously calculated that relates the GO-terms depending on their semantic similarity.
#' @param graph igraph object linking the input GO-terms.
#' @param nb_clusters Number of clusters previously calculated set to group the GO-terms.
#' @return A list containing:
#' \describe{
#' \item{graph}{An igraph object with cluster memberships assigned.}
#' \item{clusters}{A data frame with cluster assignments for each GO term.}
#' \item{nb_clusters}{The optimal number of clusters determined (if applicable).}
#' }
#' @export

clusterGO <- function(method_type = c("similarity", "network"),
                             method = NULL,
                             similarity_matrix = NULL,
                             graph, nb_clusters = NULL) {

  method_type <- match.arg(method_type)

  if (!method_type %in% c("similarity", "network")) {
    stop("method_type must be either 'similarity' or 'network'")
  }

  if (method_type == "similarity" && is.null(similarity_matrix)) {
    stop("For similarity-based clustering, a similarity_matrix is required.")
  }

  if (method_type == "similarity" && !all(V(graph)$name %in% rownames(similarity_matrix))) {
    stop("Graph vertices do not match similarity matrix names.")
  }

  if (is.null(graph)) {
    stop("For clustering, a network object is required.")
  }

  if (method_type == "similarity") {
    # Ensure similarity matrix is valid
    if (!is.matrix(similarity_matrix) || nrow(similarity_matrix) != ncol(similarity_matrix)) {
      stop("similarity_matrix must be a square matrix.")
    }

    # Perform clustering using GO semantic similarity
    dist_matrix <- as.dist(1 - similarity_matrix)

    # Perform hierarchical clustering
    hc <- hclust(dist_matrix, method = "average")

    # Cut tree to form clusters
    cluster_assignments <- cutree(hc, k = nb_clusters)

    # Assign cluster memberships to graph nodes if a graph is provided
    if (!is.null(graph)) {
      V(graph)$cluster <- cluster_assignments[V(graph)$name]
      names(cluster_assignments) <- rownames(similarity_matrix)  # Ensure correct names
    }

    return(list(graph = graph, clusters = cluster_assignments, nb_clusters = nb_clusters))

  } else if (method_type == "network") {
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop("The igraph package is required. Install it using install.packages('igraph').")
    } else {
      library(igraph)
    }
    if (!method %in% c("Walktrap", "EdgeBetweenness")) {
      stop("Invalid method for network-based clustering. Choose from 'Walktrap' or 'EdgeBetweenness'.")
    }
    if (method == "Walktrap") {
      # Perform Walktrap clustering
      clusters <- cluster_walktrap(graph)

    } else if (method == "EdgeBetweenness") {
      # Perform Edge Betweenness clustering
      clusters <- cluster_edge_betweenness(graph)
    }
    # Assign cluster memberships to nodes
    V(graph)$cluster <- clusters$membership

    return(list(graph = graph, clusters = clusters))
  }
}
