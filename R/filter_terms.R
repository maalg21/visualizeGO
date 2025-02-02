#' Filter GO-terms
#'
#' This function filters the GO-terms of the graph to obtain those that are more connected (similarity score above the threshold).
#'
#' @param similarity_matrix Squared matrix of the semantic similarity for each pair of terms calculated through any of the methods.
#' @param threshold Similarity threshold (0-1) that determines which input
#' terms will be considered not similar. Those terms for which more than half
#' do not reach the threshold will be considered not connected to the rest and
#' will be eliminated from the similarity matrix. Default = 0.1
#'
#' @return A list containing:
#' \describe{
#' \item{similarity_matrix}{A new similarity matrix without the outliers.}
#' \item{connected}{List of GO-terms that are kept in the analysis.}
#' \item{outliers}{List of GO-terms considered outliers.}
#' }
#'
#' @export

filter_terms <- function(similarity_matrix,
                         threshold = 0.1){

  # Requirements ----
  if (is.null(similarity_matrix)){
    stop("A similarity matrix is required for filtering those terms not connected (similarity score = 0).")
  }

  if (!is.matrix(similarity_matrix) ||
      nrow(similarity_matrix) != ncol(similarity_matrix)) {
    stop("Similarity matrix must be a square matrix.")
  }

  input_terms <- row.names(similarity_matrix)

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
    remove <- grep(paste(Outlier, collapse = "|"),
                   row.names(similarity_matrix))
    similarity_matrix <- similarity_matrix[-remove, -remove]
  }

  # Obtain those GO-terms with lower similarity than the threshold ----
  above_threshold <- similarity_matrix >= threshold
  above_threshold <- rowSums(above_threshold)

  result <- tryCatch({
    Connected <- which(above_threshold>=(n / 2))
  }, error = function(e) {
    stop("Error: The comparison operation failed. Lower your threshold as there is no term that has a greater similarity to it in at least half of the comparisons.")
    return(NULL)
  })

  Connected <- names(Connected)

  if (length(Connected) != 0){
    cat(length(Connected), "GO-Terms were connected above the threshold and will kept in the graph.\n")
    go_terms <- unique(Connected)
    deleted_GOterms <- input_terms[!input_terms %in% go_terms]
    cat(length(deleted_GOterms), "were removed from the graph.\n")
    remove <- grep(paste(deleted_GOterms, collapse = "|"),
                   row.names(similarity_matrix))
    similarity_matrix <- similarity_matrix[-remove, -remove]
  }

  return(list(deleted = deleted_GOterms,
              outliers = Outlier,
              connected = Connected,
              similarity_matrix = similarity_matrix))
}
