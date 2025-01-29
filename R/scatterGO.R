#' Scatter Plot of GO-terms
#'
#' Function to make a plot of GO terms as scattered points.
#' Distances between points represent the similarity between terms, and axes
#' are the first 2 components of applying a PCA to the similarity matrix.
#' Size of the point represents the number of genes the GO term contains.
#' Each color represents a cluster.
#'
#' @param similarity_matrix Squared matrix of the semantic similarity for each pair of terms calculated through any of the methods.
#' @param cluster Grouping of GO-terms in the different clusters.
#' @param OrgDb Organism to use as reference to obtain the GO-terms similarities. GOSemSimDATA object. Default = "org.Hs.eg.db" (human)
#' @param title Title of the plot
#' @param colors Vector of colors for each of the represented clusters.
#' @param labels Enable to write the Clusters number in the plot
#' @param size Determine what the determination of the size of the dots is based on.
#' It can be "genes" or "padj", if it depends on the number of genes annotated for that GO-term or on the significance value (calculated as -log10(padj)) of each GO-term respectively.
#' @param scores If the size of the GO terms represented has been chosen according to their "padj", a named numerical vector of the values to be used. Normally, this value is reflected as -log10(padj).
#'
#' @return Scatter plot of the distances of the GO-terms..
#'
#' @export
scatterGO <- function(similarity_matrix, cluster,
                      OrgDb = "org.Hs.eg.db",
                      title = NULL, colors = NULL,
                      labels = T, size = c("genes", "padj"),
                      scores = NULL) {

  if (is.null(similarity_matrix)) {
    stop("A GO semantic similarity matrix is required for performing the plot.")
  }

  if (!"clusters" %in% names(clusters) || !"graph" %in% names(clusters)) {
    stop("The 'clusters' input must contain both 'graph' and 'clusters' components.")
  }

  # Only representing the distances between inputed GO-terms
  graph <- cluster$graph
  similarity_matrix <- similarity_matrix[V(graph)$name[V(graph)$origin == "input"],
                                         V(graph)$name[V(graph)$origin == "input"]]

  distance_matrix <- 1 - similarity_matrix
  mds_result <- cmdscale(as.dist(distance_matrix), k = 2)
  pca_result <- prcomp(as.dist(distance_matrix))
  pca_scores <- pca_result$x
  pca_scores <- merge(pca_scores, by.x = "row.names",
                      data.frame("IDs" = names(cluster$clusters),
                                 "Cluster" = cluster$clusters), by.y = "IDs",
                      all.x = T) %>%
    tibble::column_to_rownames(var = "Row.names") %>%
    dplyr::mutate("Cluster" = paste("Cluster ", .$Cluster, sep = ""))

  if (size == "genes") {
    library(AnnotationDbi)
    library(OrgDb, character.only = T)

    go_terms <- c(V(graph)$name[V(graph)$origin == "input"])

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
    size <- sapply(go_gene_mapping, length)/100

  } else if (size == "padj") {
    size <- scores
  }

  library(ggplot2)
  if(is.null(colors)){
    colors <- generate_pastel_colors(n = cluster$nb_clusters)
  }

  if(isTRUE(labels)){
    library(dplyr)
    centroids <- pca_scores %>%
      group_by(Cluster) %>%
      dplyr::summarise(PC1 = mean(PC1), PC2 = mean(PC2))

    ggplot(pca_scores %>% tibble::rownames_to_column(var = "group"),
           aes(x = PC1, y = PC2, color = Cluster)) +
      geom_point(aes(size = group), alpha = 0.8) +
      geom_label(data = centroids, aes(label = Cluster),
                 vjust = -0.5, hjust = 0.5) +
      scale_color_manual(values = colors) +
      scale_size_manual(values = size) +
      geom_hline(yintercept = 0, colour = "black", linetype = "dashed") +
      geom_vline(xintercept = 0, colour = "black", linetype = "dashed") +
      theme_minimal() +
      theme(legend.position = "none",
            axis.title = element_text(size = 15, face = "bold", hjust = .5),
            plot.title = element_text(size = 20, face = "bold", hjust = .5)) +
      labs(title = title,
           x = "PC1",y = "PC2")
  } else {
    ggplot(pca_scores %>% tibble::rownames_to_column(var = "group"),
           aes(x = PC1, y = PC2, color = Cluster)) +
      geom_point(aes(size = group), alpha = 0.8) +
      scale_color_manual(values = colors) +
      scale_size_manual(values = size) +
      geom_hline(yintercept = 0, colour = "black", linetype = "dashed") +
      geom_vline(xintercept = 0, colour = "black", linetype = "dashed") +
      theme_minimal() +
      theme(legend.position = "right",
            axis.title = element_text(size = 15, face = "bold", hjust = .5),
            plot.title = element_text(size = 20, face = "bold", hjust = .5)) +
      labs(title = title,
           x = "PC1",y = "PC2")
  }
}
