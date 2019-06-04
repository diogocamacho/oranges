#' genes_pathway
genes_pathway <- function(genes_query, pathways_query) {
  gip <- apply(pathways_query, 1, function(z) paste(genes_query[which(z == 1)], collapse = " | "))
  return(gip)
}