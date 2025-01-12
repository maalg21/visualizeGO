#' Get Children
#'
#' This function get children from the GO ontology
#'
#' @param go_term List of GO-terms of interest.
#' @param go_sim_env Environment to measure semantic similarity among sets of GO terms.
#' @return All the children GO-terms of the input list.
#' @export

get_children <- function(go_term, go_sim_env) {
  # Check if the GO term exists in the children list
  if (go_term %in% names(go_sim_env$children)) {
    return(go_sim_env$children[[go_term]])
  } else {
    return(NULL)
  }
}
