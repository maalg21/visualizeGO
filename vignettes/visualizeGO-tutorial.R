## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

## ----load_data----------------------------------------------------------------
library(visualizeGO)
library(dplyr)
library(igraph)

data <- as.data.frame(readxl::read_xlsx(system.file("extdata", "Alonso-Garcia_2023.xlsx", package = "visualizeGO")))
data <- data[,c("Category", "ID", "Padj")]
GO_BP <- data[data$Category == "BP",]

## ----calculate_similarity-----------------------------------------------------
similarity_matrix <- calculate_wang(input_terms = GO_BP$ID,  ontology = "BP", # We select the category to which our GO-terms belong.
OrgDb = "org.Hs.eg.db") # We use human annotation as a reference.

## ----filter_terms-------------------------------------------------------------
filtered_terms <- filter_terms(similarity_matrix = similarity_matrix, threshold = 0.1)
similarity_matrix2 <- filtered_terms$similarity_matrix

## ----determine_nbclusters-----------------------------------------------------
set.seed(1234) # For reproducibility
# For the similarity matrix NOT filtered:
determine_nbclusters(similarity_matrix = similarity_matrix)
# Based on the Elbow method, the optimal number of clusters (k) is: 20 
# Based on the Silhouette method, the optimal number of clusters (k) is: 54

## ----determine_nbclusters2----------------------------------------------------
# For the filtered similarity matrix:
determine_nbclusters(similarity_matrix = similarity_matrix2)
# Based on the Elbow method, the optimal number of clusters (k) is: 2 
# Based on the Silhouette method, the optimal number of clusters (k) is: 34

