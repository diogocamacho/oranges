enrichment_results <- function(genes_pathway, number_genes, gene_ratio, nominal_p, p_adj) {
  enr <- tibble::tibble(data_source = pathway_info$pathway_source,
                        name = pathway_info$pathway_name,
                        number_genes = number_genes,
                        gene_ratio = gene_ratio,
                        p_val = nominal_p,
                        padj = p_adj,
                        genes = genes_pathway)
  
  enr <- enr %>% dplyr::filter(., number_genes != 0)
  
  return(enr)
}