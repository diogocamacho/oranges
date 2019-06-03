#' Plot enrichment results
#'
#' \code{plot_oranges} returns a plot for the top 25 enriched pathways identified with the \code{\link{oranges}} package.
#'
#' @param res ORANGES results data frame
#' @param thr Significance (FDR-corrected q-value) for pathway enrichment.  Defaults to q < 0.05
plot_oranges <- function(res, num_pathways) {
  
  if(missing(num_pathways)) num_pathways <- 100
  
  if(length(unique(res$data_source)) > 1) {
    res %>% 
      dplyr::filter(., padj != 1) %>%
      dplyr::arrange(., padj) %>%
      dplyr::slice(1:num_pathways) %>%
      ggplot() + 
      geom_point(aes(x = name, y = -log10(padj), color = gene_overlap, size = number_genes), alpha = 0.5) +
      scale_color_viridis_c() + 
      facet_grid(. ~ data_source, scales = "free") + 
      theme_bw() + 
      theme(axis.text.x = element_blank(),
            panel.grid = element_blank())
  } else {
    res %>% 
      dplyr::filter(., padj != 1) %>%
      dplyr::arrange(., padj) %>%
      dplyr::slice(1:num_pathways) %>%
      ggplot() + 
      geom_point(aes(x = name, y = -log10(padj), color = gene_overlap, size = number_genes), alpha = 0.5) +
      scale_color_viridis_c() + 
      theme_bw() + 
      theme(axis.text.x = element_blank(),
            panel.grid = element_blank())
  }
}

