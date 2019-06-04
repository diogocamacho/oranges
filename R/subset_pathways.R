subset_pathways <- function(query_entrez, universe_entrez) {
  x1 <- which(colnames(cpdb_data$pathway_matrix) %in% universe_entrez)
  x2 <- which(colnames(cpdb_data$pathway_matrix) %in% query_entrez)
  ids <- list(uid = x1, qid = x2)
  return(ids)
}