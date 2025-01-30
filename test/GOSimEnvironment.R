# Load necessary libraries
library(AnnotationDbi)
library(GO.db)

# Initialize the database list
GOSimEnv <- new.env()

# Define the relationships for children and parents for all ontologies
children_function <- list(
  "BP" = GOBPCHILDREN,
  "CC" = GOCCCHILDREN,
  "MF" = GOMFCHILDREN
)

parents_function <- list(
  "BP" = GOBPPARENTS,
  "CC" = GOCCPARENTS,
  "MF" = GOMFPARENTS
)

# Retrieve all GO terms for all ontologies (BP, CC, MF)
all_go_terms <- unique(c(keys(GO.db, keytype = "GOID", ontology = "BP"),
                         keys(GO.db, keytype = "GOID", ontology = "CC"),
                         keys(GO.db, keytype = "GOID", ontology = "MF")))

# Initialize lists for ancestors and children
ancestors_list <- list()
children_list <- list()

# Iterate over all GO terms (assuming 'all_go_terms' is the list of GO terms you want to process)
for (go_term in all_go_terms) {
  
  # Get the ontology of the current GO term using the GO.db annotation data
  ont_info <- AnnotationDbi::select(GO.db, keys = go_term, columns = "ONTOLOGY", keytype = "GOID")
  
  # Check if the ontology information is available for the GO term
  if (nrow(ont_info) > 0) {
    ont <- ont_info$ONTOLOGY[1]  # Get the ontology (BP, CC, or MF)
    
    # Use the appropriate function for children and parents based on the ontology
    if (ont == "BP") {
      children <- GOBPCHILDREN[[go_term]]
      ancestors <- GOBPPARENTS[[go_term]]
    } else if (ont == "CC") {
      children <- GOCCCHILDREN[[go_term]]
      ancestors <- GOCCPARENTS[[go_term]]
    } else if (ont == "MF") {
      children <- GOMFCHILDREN[[go_term]]
      ancestors <- GOMFPARENTS[[go_term]]
    }
    
    # Store the children and ancestors in the respective lists
    children_list[[go_term]] <- children
    ancestors_list[[go_term]] <- ancestors
  }
}

# Store the lists in the environment as $ancestors and $children
GOSimEnv$ancestors <- ancestors_list
GOSimEnv$children <- children_list

# Save the GOSimEnv to an .RData file
save(GOSimEnv, file = "go_term_database_all_ontologies.RData")

# The go_term_database_all_ontologies.RData file is now created and contains the ancestors and children lists for all GO terms.