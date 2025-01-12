#' Transform GO IDs into numbers
#'
#' This function transform GO IDs (GO:00000000) into numbers to make the plot more visible.
#'
#' @param graph The graph that is going to be plotted.
#' @return Return a data frame of the corresponding number and GO ID.
#' @export

transform_go_ids_to_numbers <- function(graph) {
  # Create a mapping of GO term IDs to numbers
  go_terms <- V(graph)$name
  term_to_number <- setNames(seq_along(go_terms), go_terms)

  # Update vertex names to numeric IDs
  V(graph)$name <- as.character(sapply(V(graph)$name, function(go) term_to_number[go]))

  # Create a data frame mapping GO ID to number and description (if available)
  # Fetch descriptions for GO terms from GO.db
  go_descriptions <- sapply(go_terms, function(go_term) {
    # Fetch description using GO.db
    desc <- tryCatch({
      Term(go_term)  # Get description using Term function from GO.db
    }, error = function(e) {
      NA  # If not found, return NA
    })
    return(desc)
  })

  # Create a data frame mapping GO ID to number and description
  GODescriptions <- data.frame(
    GO_ID = names(term_to_number),
    Number = term_to_number,
    Description = go_descriptions,
    stringsAsFactors = FALSE
  )

  return(list(graph = graph, GODescriptions = GODescriptions))
}
