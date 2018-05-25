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
signature_enrichment <- function(query_set,
                      query_universe,
                      query_gene_names,
                      signature_geneset,
                      signature_universe,
                      signature_source,
                      signature_name)
{
  
  # Cleaning NA in query set
  if(length(which(is.na(query_set))) != 0)
  {
    query_set <- query_set[-which(is.na(query_set))]
  }
  
  # define signature
  signature_enrichment <- matrix(0,1,5)
  #signature_enrichment <- data.frame()
  nix <- matrix(0,1,1)
  my_universe <- intersect(unique(query_universe),signature_universe) # elements of genes where query set came from that are in the universe to be enriched
  M <- length(my_universe) # elements of genes where query set came from that are in the universe to be enriched
  
  # compute p-value
  clean.set <- intersect(signature_geneset[which(!is.na(signature_geneset))],my_universe)
  absent.in.set <- setdiff(query_set,clean.set)
  refined.query <- setdiff(query_set,absent.in.set)
  N <- length(refined.query)
  X <- length(intersect(refined.query,clean.set))
  X.genes <- paste(as.character(query_gene_names[match(intersect(refined.query,clean.set),query_universe)]),collapse=" | ")
  K <- length(clean.set)
  a1 <- X
  a2 <- N - a1
  a3 <- K - a1 
  a4 <- M - (a1 + a2 + a3)
  
  if (N == 0)
  {
    signature_enrichment <- cbind(as.character(signature_source),as.character(signature_name),0,1,"NA")
    # signature_enrichment$signature_source <- signature_source
    # signature_enrichment$signature_name <- signature_name
    # signature_enrichment$representation <- paste(N,"/",K,sep="")
    # signature_enrichment$nominal_pvalue <- 1
    # signature_enrichment$genes_in_pathway <- NA
  } else {
    tmp.res <- fisher.test(cbind(c(a1,a3),c(a2,a4)),alternative="two.sided")$p.value
    signature_enrichment <- cbind(as.character(signature_source),as.character(signature_name),N/K,tmp.res,X.genes)
    # signature_enrichment$signature_source <- signature_source
    # signature_enrichment$signature_name <- signature_name
    # signature_enrichment$representation <- paste(N,"/",K,sep="")
    # signature_enrichment$nominal_pvalue <- tmp.res
    # signature_enrichment$genes_in_pathway <- X.genes
  }  
  
  return(signature_enrichment)
  
}