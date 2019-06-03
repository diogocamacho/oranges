extract_pathways <- function(ids) {
  unames <- gene_names[ids$uid]
  qnames <- gene_names[ids$qid]
  
  sub_p <- P[, ids$uid]
  sub_q <- P[, ids$qid]
  
  eps <- list(genes_universe = unames,
              genes_query = qnames,
              pathways_universe = sub_p,
              pathways_query = sub_q)
  
  return(eps)
}
