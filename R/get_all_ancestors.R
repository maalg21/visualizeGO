#' Get Ancestors
#'
#' This function get all the ancestors from the input GO-terms list
#'
#' @param go_list List of GO-terms of interest.
#' @param go_sim_env Environment to measure semantic similarity among sets of GO terms.
#' @return All the ancestors GO-terms of the input list.
#' @export

get_all_ancestors <- function(go_list, go_sim_env) {
  all_ancestors <- unique(unlist(lapply(go_list, function(go_term) {
    # Check if the GO term exists in the ancestors list
    if (go_term %in% names(go_sim_env$ancestors)) {
      return(go_sim_env$ancestors[[go_term]])
    } else {
      return(NULL)
    }
  })))

  return(all_ancestors)
}
