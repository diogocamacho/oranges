#' Fisher's exact test for pathway enrichment
#'
#' \code{ora} returns the enrichment of pathway sets given a query gene signature based on Fisher's exact test. Main function in the \code{\link{oranges}} package.
#'
#' @param query_set Gene signature of interest, as EntrezIDs
#' @param universe_entrez Full EntrezID set where the query set comes from (usually the full set of RNA-seq or microarray genes characterized in the transcriptomics platform)
#' @param universe_symbols Gene symbols for the entire set of genes in universe
#' @param pathway_matrix A (sparse) matrix where each row is a pathway and each column a gene. As distributed, this matrix has pathway sets available in CPDB and MSigDB.
#' @return A tibble with overlap of query set with pathway genes, set of query genes in pathway (when found), and Fisher's exact test p-value
ora <- function(query_set, universe_entrez, universe_symbols, pathway_matrix)
{
  
  # Cleaning NA in query set
  if (length(which(is.na(query_set))) != 0) query_set <- query_set[-which(is.na(query_set))]
  
  # fisher's exact test for each signature
  enr <- vector(mode = "list",length = nrow(pathway_matrix))
  
  N <- length(query_set) # <-- how big query set is
  K <- Matrix::rowSums(pathway_matrix) # <-- sizes of pathways
  M <- ncol(pathway_matrix) # <-- size of my universe
  
  for (i in seq(1,length(enr))) {
    
    # all genes in pathway
    x1 <- names(which(pathway_matrix[i, ] == 1)) 
    # M <- ncol(pathway_matrix) # <-- size of my universe
    
    # confusion matrix
    # N <- length(query_set) # <-- how big query set is
    
    x2 <- intersect(query_set, x1)
    X <- length(x2)
    X_symbols <- paste(as.character(universe_symbols[which(universe_entrez %in% x2)]), collapse=" | ")
    
    # K <- length(x1) # <-- number of genes in pathway
    
    a2 <- N - X
    a3 <- K[i] - X 
    a4 <- M - (X + a2 + a3)
    
    # compute fisher's exact test
    if (N == 0)
    {
      enr[[i]] <- tibble::tibble(gene_overlap = 0,
                                 genes_pathway = NA,
                                 p_val = 1)
    } else {
      res <- fisher.test(cbind(c(X, a3), c(a2, a4)), alternative = "two.sided")$p.value
      enr[[i]] <- tibble::tibble(gene_overlap = X/K[i],
                                 genes_pathway = X_symbols,
                                 p_val = res)
    }
  }
  enr <- dplyr::bind_rows(enr)
  return(enr)
}