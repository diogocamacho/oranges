extract_pathways <- function(ids) {
  unames <- cpdb_data$gene_symbols[ids$uid]
  qnames <- cpdb_data$gene_symbols[ids$qid]
  
  sub_p <- cpdb_data$pathway_matrix[, ids$uid]
  sub_q <- cpdb_data$pathway_matrix[, ids$qid]
  
  eps <- list(genes_universe = unames,
              genes_query = qnames,
              pathways_universe = sub_p,
              pathways_query = sub_q)
  
  return(eps)
}
