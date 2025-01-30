#' Filter GO-terms
#'
#' This function filters the GO-terms of the graph to obtain those that are more connected (similarity score above the threshold).
#'
#' @param similarity_matrix Squared matrix of the semantic similarity for each pair of terms calculated through any of the methods.
#' @param cluster Grouping of GO-terms in the different clusters.
#' @param threshold Similarity threshold (0-1). Default = 0.1
#' @param input_terms Important GO-terms (tipically those that you analyzed) that you want to keep in the next steps. Normally, only those introduced at the start are kept.
#'
#' @return A list containing:
#' \describe{
#' \item{graph}{A filtered igraph object, where the least connected terms have been eliminated.}
#' \item{similarity_matrix}{}
#' \item{deleted}{List of GO-terms eliminated for not presenting a similarity greater than the threshold in at least half of the relationships.}
#' \item{connected}{List of GO-terms that are kept in the analysis.}
#' \item{outliers}{List of GO-terms considered outliers.}
#' }
#'
#' @export

filter_terms <- function(similarity_matrix,
                         graph, threshold = 0.1,
                         input_terms = NULL){

  # Requirements ----
  if (is.null(graph)) {
    stop("For clustering, a network object is required.")
  }

  if (is.null(similarity_matrix)){
    stop("A similarity matrix is required for filtering those terms not connected (similarity score = 0).")
  }

  if (!is.matrix(similarity_matrix) || nrow(similarity_matrix) != ncol(similarity_matrix)) {
    stop("Similarity matrix must be a square matrix.")
  }

  # Obtain those terms with similarity score = 0 ----
  # Count zeros in each row
  zero_counts <- rowSums(similarity_matrix == 0)

  # Matrix size
  n <- nrow(similarity_matrix)

  # Get row indices where at least half of the values are zero
  Outlier <- which(zero_counts == n - 1)
  Outlier <- names(Outlier)

  if (length(Outlier) != 0){
    cat("GO-Terms", Outlier, "were not similar to any other GO-term. They are outliers.\n")
    cat("It'll be removed from the graph.\n")
    graph <- delete_vertices(graph, Outlier)
    remove <- grep(Outlier, row.names(similarity_matrix))
    similarity_matrix <- similarity_matrix[-remove, -remove]
  }

  # Obtain those GO-terms with lower similarity than the threshold ----
  above_threshold <- similarity_matrix >= threshold
  above_threshold <- rowSums(above_threshold)
  Connected <- which(above_threshold>=(n / 2))
  Conected <- names(Connected)

  if (length(Connected) != 0){
    cat(length(Connected), "GO-Terms were connected above the threshold and will kept in the graph.\n")
    cat(length(Conected[!V(graph)[V(graph)$origin == "input"]$name %in% Conected]), "are input GO-terms that were not connected above the threshold BUT those terms will be kept.\n")
    go_terms <- unique(c(V(graph)[V(graph)$origin == "input"]$name, Conected))
    deleted_GOterms <- V(graph)[!V(graph)$name %in% go_terms]$name
    graph <- delete_vertices(graph, deleted_GOterms)
    cat(paste(deleted_GOterms, collapse = ", "), "were removed from the graph.\n")
    remove <- grep(paste(deleted_GOterms, collapse = "|"),
                   row.names(similarity_matrix))
    similarity_matrix <- similarity_matrix[-remove, -remove]
  }

  return(list(graph = graph,
              deleted = deleted_GOterms,
              outliers = Outlier,
              connected = Connected,
              similarity_matrix = similarity_matrix))
}
