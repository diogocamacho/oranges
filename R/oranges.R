#' ORANGES: Over-Representation ANalysis for Gene Expression Signatures
#'
#'
#' \code{oranges} returns the enrichment of pathway sets given a query gene signature based on Fisher's exact test. Currently this function is built around enrichment of human pathway sets, but \code{\link{ora}} function can be generalized to fit different organisms of interest provided a pathway matrix is provided.
#'
#' @param query_set Gene signature of interest, as EntrezIDs
#' @param universe_entrez Full EntrezID set where the query set comes from (usually the full set of RNA-seq or microarray genes characterized in the transcriptomics platform)
#' @param universe_symbols Gene symbols for the entire set of genes in universe
#' @return A tibble.
oranges <- function(query_entrez, universe_entrez) {

  wp <- subset_pathways(query_entrez = query_entrez, 
                        universe_entrez = universe_entrez) #<-- returns ids (columns)

  eps <- extract_pathways(wp)

  gid <- genes_pathway(genes_query = eps$genes_query, 
                       pathways_query = eps$pathways_query)

  enr <- calculate_enrichment(genes_query = wp$genes_query,
                              pathways_universe = wp$pathways_universe,
                              pathways_query = wp$pathways_query)

  res <- enrichment_results(genes_pathway = gid,
                            number_genes = enr$number_genes,
                            gene_ratio = enr$gene_ratio,
                            nominal_p = enr$nominal_p,
                            p_adj = enr$padj)

  return(res)
}
