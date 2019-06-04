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
      geom_point(aes(x = name, y = -log10(padj), color = pathway_proportion, size = number_genes), alpha = 0.5) +
      scale_color_viridis_c() + 
      labs(x = "Pathway set", y = "-log10(adjusted p-value)") +
      facet_grid(. ~ data_source, scales = "free") + 
      theme_bw() + 
      theme(axis.ticks.length = unit(0.5, "cm"),
            axis.ticks = element_line(size = 1),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 12, color = "black"),
            axis.title = element_text(size = 12, color = "black"),
            panel.grid = element_blank(),
            strip.text.x = element_text(size = 12, color = "black"))
  } else {
    res %>% 
      dplyr::filter(., padj != 1) %>%
      dplyr::arrange(., padj) %>%
      dplyr::slice(1:num_pathways) %>%
      ggplot() + 
      geom_point(aes(x = name, y = -log10(padj), color = pathway_proportion, size = number_genes), alpha = 0.5) +
      scale_color_viridis_c() + 
      labs(x = "Pathway set", y = "-log10(adjusted p-value)") +
      theme_bw() + 
      theme(axis.ticks.length = unit(0.5, "cm"),
            axis.ticks = element_line(size = 1),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 12, color = "black"),
            axis.title = element_text(size = 12, color = "black"),
            panel.grid = element_blank(),
            strip.text.x = element_text(size = 12, color = "black"))
  }
}

