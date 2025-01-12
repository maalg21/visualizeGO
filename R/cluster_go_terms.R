#' Cluster GO Terms
#'
#' This function performs clustering of Gene Ontology (GO) terms based on either semantic similarity or network-based methods.
#'
#' @param method_type Choose between "similarity" or "network" to select the method to assign the clusters to each input GO term.
#' @param method Once selected the type of assignment, it's the method within each type selected to obtain the clusters.
#' @param ontology Gene Ontology category to use (could be "BP" for Biological Process, "CC" for Cellular Component or "MF" for "Molecular Function").
#' @param orgdb Organism to use as reference to obtain the GO-terms similarities. Default: "org.Hs.eg.db"
#' @param similarity_matrix Only if the method_type selected was "similarity". Similarity matrix previously calculated that relates the GO-terms depending on their semantic similarity.
#' @param graph igraph object linking the input GO-terms.
#' @param nb_clusters If it's known, number of clusters set to group the GO-terms.
#' @param k_range Specified range to uses silhouette scores to determine the optimal number of clusters from.
#' @return A list containing:
#' \describe{
#' \item{graph}{An igraph object with cluster memberships assigned.}
#' \item{clusters}{A data frame with cluster assignments for each GO term.}
#' \item{nb_clusters}{The optimal number of clusters determined (if applicable).}
#' }
#' @export

cluster_go_terms <- function(method_type = c("similarity", "network"),
                             method = NULL,
                             orgdb = "org.Hs.eg.db",
                             ontology = c("BP", "CC", "MF"),
                             similarity_matrix = NULL,
                             graph,
                             nb_clusters = NULL,
                             k_range = 2:10) {

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

    # Assign cluster memberships to graph nodes
    if(is.null(nb_clusters)){
      # Calculate silhouette scores for a range of k
      silhouette_scores <- sapply(k_range, function(k) {
        cluster_assignment <- cutree(hc, k = k)
        sil <- cluster::silhouette(cluster_assignment, dist_matrix)
        mean(sil[, 3])  # Extract the average silhouette width
      })

      # Find the optimal k based on maximum silhouette score
      optimal_k <- k_range[which.max(silhouette_scores)]
      print(paste("Optimal k using Silhouette Method:", optimal_k))
      nb_clusters <- optimal_k
    }

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
      stop("The igraph package is required. Install it using install.packages('GOSemSim').")
    } else {
      library(igraph)
    }
    if (!method %in% c("Walktrap", "Louvain", "EdgeBetweenness")) {
      stop("Invalid method for network-based clustering. Choose from 'Walktrap', 'Louvain', or 'EdgeBetweenness'.")
    }
    if (method == "Walktrap") {
      # Perform Walktrap clustering
      clusters <- cluster_walktrap(graph)

    } else if (method == "Louvain") {
      # Perform Louvain clustering
      clusters <- cluster_louvain(graph)

    } else if (method == "EdgeBetweenness") {
      # Perform Edge Betweenness clustering
      clusters <- cluster_edge_betweenness(graph)
    }
    # Assign cluster memberships to nodes
    V(graph)$cluster <- clusters$membership

    return(list(graph = graph, clusters = clusters))
  }
}
