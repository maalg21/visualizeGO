#' @title Build Hierarchical Graph
#'
#' @description This function links all parents and children in the list of input GO terms and creates a graph, establishing hierarchical relationships between the GO terms.
#'
#' @param cluster Grouping of GO-terms in the different clusters.
#' @param size Determine what the determination of the size of the dots is based on.
#' It can be "genes" or "padj", if it depends on the number of genes annotated for that GO-term or on the significance value (calculated as -log10(padj)) of each GO-term respectively.
#' @param OrgDb Organism to use as reference to obtain the GO-terms similarities. GOSemSimDATA object. Default = "org.Hs.eg.db" (human)
#' @param scores If the size of the GO terms represented has been chosen according to their "padj", a named numerical vector of the values to be used. Tipically, this value is reflected as -log10(padj).
#' @param title Title of the plot
#' @param colors Vector of colors for each of the represented clusters.
#' @return A tree map plot relating all the input GO terms to their clusters.
#' @export

treeMap <- function(cluster, size, scores,
                    OrgDb = "org.Hs.eg.db",
                    colors = NULL, title = NULL){

  # Ensure all needed packages are installed
  if (!requireNamespace("treemap", quietly = TRUE)) {
    stop("The 'treemap' package is required but not installed. Please install it using install.packages('treemap').")
  }
  library(treemap)

  if (!requireNamespace("treemapify", quietly = TRUE)) {
    stop("The 'treemapify' package is required but not installed. Please install it using install.packages('treemapify').")
  }
  library(treemapify)

  library(ggplot2)

  go_data <- merge(cluster$representative_pathways,
                   cluster$descriptions %>%
                     mutate("Cluster" = paste("Cluster ", Cluster, sep = "")),
                   by = "Cluster")

  if(size == "genes"){
    cat("This option takes way long time to process ... Be patient!")
    library(AnnotationDbi)
    library(OrgDb, character.only = T)

    go_terms <- go_data$GO_ID

    # Get genes annotated to the GO terms
    go_gene_mapping <- list()
    for(g in 1:length(go_terms)){
      tmp <- AnnotationDbi::select(x = get(OrgDb),
                                   keys = go_terms[g],
                                   keytype = "GOALL",
                                   columns = "SYMBOL")
      genes <- c(tmp$SYMBOL)
      go_gene_mapping[[g]] <- unique(genes)
    }
    names(go_gene_mapping) <- go_terms
    scores <- sapply(go_gene_mapping, length)/100
  }

  go_data <- merge(go_data, as.data.frame(scores) %>%
                       tibble::rownames_to_column(var = "GO_ID"),
                   by = "GO_ID") %>%
      dplyr::select(-GO.IDs)

  if(is.null(colors)){
    colors <- generate_pastel_colors(n = length(unique(cluster$clusters)))
  }

  ggplot(go_data, aes(area = scores,
                      fill = Cluster, label = Description,
                      subgroup = Representative.Pathway)) +
    geom_treemap(alpha = .5) +
    geom_treemap_text(colour = "white",
                      size = 8, alpha = 0.7,
                      place = "center") +
    geom_treemap_subgroup_text(place = "center",
                               size = 12, colour = "black",
                               fontface = "bold", alpha = 1) +
    scale_fill_manual(values = colors) +
    labs(title = title) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold",
                                    hjust = .5, size = 15),
          plot.tag = element_text(face = "bold", size = 12))

}
