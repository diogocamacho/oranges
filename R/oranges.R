#' ORANGES: Over-Representation ANalysis for Gene Expression Signatures
#' 
#'
#' \code{oranges} returns the enrichment of pathway sets given a query gene signature based on Fisher's exact test. Currently this function is built around enrichment of human pathway sets, but \code{\link{ora}} function can be generalized to fit different organisms of interest provided a pathway matrix is provided.
#'
#' @param query_set Gene signature of interest, as EntrezIDs
#' @param universe_entrez Full EntrezID set where the query set comes from (usually the full set of RNA-seq or microarray genes characterized in the transcriptomics platform)
#' @param universe_symbols Gene symbols for the entire set of genes in universe
#' @return A tibble.
oranges <- function(query_entrez,universe_entrez,universe_names)
{
  # load pathway sets
  load("./data/pathway_sets.RData")

  # pathway matrix
  P <- pathway_sets[[1]]
  
  # pathway metadata
  pathway_info <- pathway_sets[[2]]
  
  # clean up matrix & universe genes
  keep_cols <- which(colnames(P) %in% universe_entrez)
  keep_ids <- which(universe_entrez %in% colnames(P)[keep_cols])
  universe_entrez <- universe_entrez[keep_ids]
  universe_symbols <- universe_symbols[keep_ids]
  P <- P[,keep_cols]
  
  # clean matrix
  nix <- Matrix::rowSums(P)
  nix <- which(nix < 3 & nix > 500)
  if (length(nix) != 0) {
    P <- P[-nix, ]
    pathway_info <- pathway_info[-nix, ]
  }
  
  # perform enrichment
  enr <- ora(query_set = query_entrez,
             universe_entrez = universe_entrez,
             universe_symbols = universe_symbols,pathway_matrix = P)
  
  enr <- enr %>% 
    add_column(.,pathway_source = pathway_info$pathway_source,.before = 1) %>% 
    add_column(.,pathway_name = pathway_info$pathway_name,.before = 2) %>%
    add_column(.,q_val = p.adjust(enr$p_val,"fdr")) %>% 
    dplyr::filter(.,gene_overlap != 0)

  return(enr)
}