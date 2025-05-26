#' Generate Clusters Table
#'
#' The generate_cluster_table function creates a detailed and visually appealing table
#' summarizing the clustering results of Gene Ontology (GO) terms. The function can work
#' with both similarity-based clustering and network-based clustering outputs.
#'
#' This function is useful for reporting clustering results in a clear, interactive,
#' and visually appealing format, ideal for presentations or publications.
#'
#' @description
#' A formatted HTML table summarizing: (1) Cluster ID, (2) Representative Pathway
#' (for similarity-based clustering) and (3) GO IDs within each cluster.
#' Adds row-specific colors to highlight clusters visually.
#'
#' @param cluster_output Result from "clusterGO".
#' @param col_palette Colors to use for each cluster.
#' @param text_color Color of the text in the table.
#' @param file_name Name of the saved file.
#' @return Combines cluster data into a single data frame.
#' Assigns colors to rows based on cluster memberships using the provided or
#' default color palette. If enabled, saves the table as a PNG file.
#' @export

generate_cluster_table <- function(cluster_output,
                                   col_palette = NULL,
                                   text_color = "black",
                                   file_name = NULL,
                                   width = 800,
                                   height = 600,
                                   zoom = 2) {

  if (is.null(cluster_output) ||
      !"clusters" %in% names(cluster_output) ||
      !"descriptions" %in% names(cluster_output) ||
      !"representative_pathways" %in% names(cluster_output)) {
    stop("Invalid cluster result provided.")
  }

  cluster_df <- cluster_output$representative_pathways

  # Generate a table image
  # Use kableExtra for better table formatting
  if (!requireNamespace("kableExtra", quietly = TRUE)) {
    stop("The kableExtra package is required for table visualization. Install it using install.packages('kableExtra').")
  }
  library(kableExtra)

  # Add colors to rows based on clusters
  # Visualize clusters (assign colors)
  if(is.null(col_palette)){
    cluster_colors <- generate_pastel_colors(n = length(unique(cluster_output$clusters)))
  } else {
    cluster_colors <- col_palette
  }

  names(cluster_colors) <- paste("Cluster ", unique(cluster_output$clusters), sep = "")

  final_table <- cluster_df %>%
    kbl(format = "html", escape = FALSE, row.names = FALSE) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = FALSE) %>%
    row_spec(0, bold = TRUE) %>%
    column_spec(1, color = text_color,
                background = cluster_colors[1:nrow(cluster_df)])

  print(final_table)

  if(!is.null(file_name)){
    temp_html <- paste(file_name, ".html")

    cluster_df %>%
      kbl(format = "html", escape = FALSE, row.names = FALSE) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                    full_width = FALSE) %>%
      row_spec(0, bold = TRUE) %>%
      column_spec(1, color = text_color,
                  background = cluster_colors[1:nrow(cluster_df)]) %>%
      save_kable(temp_html)  # Save the table to an HTML file
    # Use webshot to convert HTML to PNG
    library(webshot)
    webshot(temp_html, file_name, vwidth = width, vheight = height, zoom = zoom)
  }

  invisible(final_table)

}
