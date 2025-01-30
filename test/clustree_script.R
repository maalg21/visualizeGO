dist_matrix <- as.dist(1 - similarity_matrix)
hc <- hclust(dist_matrix, method = "average")

tmp <- list()
for (k in 1:55){
  tmp[[k]] <- cutree(hc, k = k)
}
df <- as.data.frame(do.call(cbind, tmp))
colnames(df) <- paste("K", 1:55, sep = "")

pca <- prcomp(dist_matrix, center = TRUE, scale. = FALSE)
ind.coord <- pca$x

ind.coord <- merge(ind.coord, by.x = "row.names",
                   as.data.frame(df) %>%
                     tibble::rownames_to_column(var = "GOTerm"),
                   by.y = "GOTerm", all = T) %>%
  tibble::column_to_rownames(var = "Row.names")
clustree::clustree(ind.coord, prefix = "K")
