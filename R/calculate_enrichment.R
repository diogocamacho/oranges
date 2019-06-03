calculate_enrichment <- function(genes_query, pathways_universe, pathways_query) {
  N <- q_length # <-- how big query set is
  K <- Matrix::rowSums(pathways_universe) # <-- sizes of pathways
  M <- ncol(pathways_universe) # <-- size of my universe
  
  # X <- length(x2) # ----> subq
  X <- Matrix::rowSums(pathways_query)
  
  a2 <- N - X
  a3 <- K - X 
  a4 <- M - (X + a2 + a3)
  
  FT <- cbind(X = X, fn = a3, fp = a2, tn = a4)
  sig_enr <- apply(FT, 1, function(x) fisher.test(matrix(x, 2, 2), alternative = "two.sided")$p.value)
  
  ug <- unique(pathway_info$group)
  cor_p <- integer(length = length(sig_enr))
  for (i in seq(1, length(ug))) {
    cor_p[pathway_info$group == ug[i]] <- p.adjust(sig_enr[pathway_info$group == ug[i]], "fdr")
  }
  
  stats_enr <- list(number_genes = X,
                    gene_ratio = paste(as.character(X), "/", as.character(K), sep = ""),
                    nominal_p = sig_enr,
                    padj = cor_p)
  return(stats_enr)
}