#' Cluster GO Terms
#'
#' This function performs clustering of Gene Ontology (GO) terms based on either semantic similarity or network-based methods.
#'
#' @param similarity_matrix Similarity matrix previously calculated that relates the GO-terms depending on their semantic similarity.
#' @param nb_clusters Number of clusters previously calculated set to group the GO-terms.
#' @return A list containing:
#' \describe{
#' \item{clusters}{A data frame with cluster assignments for each GO term.}
#' }
#' @export

clusterGO <- function(similarity_matrix = NULL,
                      nb_clusters = NULL) {

  if (is.null(similarity_matrix)) {
    stop("For similarity-based clustering, a similarity_matrix is required.")
  }

  if (!is.matrix(similarity_matrix) || nrow(similarity_matrix) != ncol(similarity_matrix)) {
      stop("similarity_matrix must be a square matrix.")
  }

  # Save input terms ID
  input_terms <- row.names(similarity_matrix)

  # Perform clustering using GO semantic similarity
  dist_matrix <- as.dist(1 - similarity_matrix)

  # Perform hierarchical clustering
  hc <- hclust(dist_matrix, method = "average")

  # Cut tree to form clusters
  cluster_assignments <- cutree(hc, k = nb_clusters)

  # Helper function to get GO term names from GO IDs
  get_go_term_name <- function(go_id) {
    go_term_name <- tryCatch({
      AnnotationDbi::Term(go_id)
    }, error = function(e) {
      return(go_id)  # If the GO term is not found, return the GO ID itself
    })
    return(go_term_name)
  }

  go_descriptions <- c()
  # Obtain the term name
  for(g in input_terms){
    tmp <- get_go_term_name(g)
    if(is.na(tmp)){
      names(tmp) <- g
    }
    go_descriptions <- c(go_descriptions, tmp)
  }

  # Create a data frame mapping GO ID to number and description
  GODescriptions <- data.frame(
    GO_ID = input_terms,
    Description = go_descriptions,
    stringsAsFactors = FALSE
  )

  GODescriptions <- merge(GODescriptions, by = "GO_ID",
                          as.data.frame(cluster_assignments) %>%
                            tibble::rownames_to_column(var = "GO_ID") %>%
                            dplyr::rename("Cluster" = cluster_assignments)) %>%
    arrange(Cluster)
  # This stores the GO Descriptions data frame

  # For each cluster, determine the most representative term and list GO IDs
  clusters <- unique(GODescriptions$Cluster)

  # Initialize a data frame for the table
  cluster_info <- list()

  for (cluster_id in unique(clusters)) {
    go_ids <- names(which(cluster_assignments == cluster_id))

    # Calculate the most representative term
    # Representative term: the term with the highest average similarity to other terms in the cluster
    sub_matrix <- similarity_matrix[go_ids, go_ids, drop = FALSE] # Subset similarity matrix
    avg_similarity <- rowMeans(sub_matrix)
    representative_term_id <- names(which.max(avg_similarity))

    # Get the full name of the representative term
    representative_term_name <- get_go_term_name(representative_term_id)

    # Check for NA values before proceeding
    if (is.na(representative_term_name) || length(go_ids) == 0) {
      cat("In Cluster", cluster_id, ", no representative patwhay was found.")
      representative_term_name <- representative_term_id
    }

    # Store information in the list
    cluster_info[[paste("Cluster", cluster_id)]] <- data.frame(
      Cluster = paste("Cluster", cluster_id),
      "Representative Term" = representative_term_name,
      "GO IDs" = paste(go_ids, collapse = ", "),
      stringsAsFactors = FALSE
    )
  }

  # Combine all cluster info into a single data frame
  cluster_df <- do.call(rbind, cluster_info)

  return(list(clusters = cluster_assignments,
              descriptions = GODescriptions,
              representative_term_name = cluster_df))
}
