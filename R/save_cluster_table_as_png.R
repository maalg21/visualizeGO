#' Save Clusters Table as PNG
#'
#'
#' @param cluster_df Output of "generate_cluster_table" function.
#' @param file_name Name of the saved file.
#' @return Saves the table as a PNG image.
#' @seealso [generate_cluster_table()]
#' @export

save_cluster_table_as_png <- function(cluster_df, file_name = "cluster_table.png",
                                      width = 800, height = 600, zoom = 2) {

  # Save the table as an HTML file
  temp_html <- tempfile(fileext = ".html")

  # Generate the HTML table using kableExtra
  cluster_df %>%
    kbl(format = "html", escape = FALSE, row.names = FALSE) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = FALSE) %>%
    save_kable(temp_html)  # Save the table to an HTML file

  # Use webshot to convert HTML to PNG
  library(webshot)
  webshot(temp_html, file_name, vwidth = width, vheight = height, zoom = zoom)
}
