##############################################################################
# SIGNATURE ENRICHMENT
# This will compute the p-value based on Fisher's test to determine enrichment 
# of a gene set against a gene signature.  Takes in a set of entrez ids.
#
# be aware of the universes: the query set universe is the one where the query set
# was derived from, which is not necessarily the same as the signature universe
#
# signature.source and signature.name are character entries
##############################################################################
ora <- function(query_set,universe_entrez,universe_symbols,pathway_matrix)
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
    x1 <- names(which(pathway_matrix[i,] == 1)) 
    # M <- ncol(pathway_matrix) # <-- size of my universe
    
    # confusion matrix
    # N <- length(query_set) # <-- how big query set is
    
    x2 <- intersect(query_set,x1)
    X <- length(x2)
    X_symbols <- paste(as.character(universe_symbols[which(universe_entrez %in% x2)]),collapse=" | ")
    
    # K <- length(x1) # <-- number of genes in pathway
    
    a1 <- X
    a2 <- N - a1
    a3 <- K[i] - a1 
    a4 <- M - (a1 + a2 + a3)
    
    # compute fisher's exact test
    if (N == 0)
    {
      enr[[i]] <- data_frame(gene_overlap = 0,genes_pathway = NA,p_val = 1)
    } else {
      res <- fisher.test(cbind(c(a1,a3),c(a2,a4)),alternative="two.sided")$p.value
      enr[[i]] <- data_frame(gene_overlap = X/K[i],genes_pathway = X_symbols,p_val = res)
    }
  }
  enr <- bind_rows(enr)
  return(enr)
}