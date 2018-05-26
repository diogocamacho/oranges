#' Plot enrichment results
#'
#' \code{plot_oranges} returns a plot for the top 25 enriched pathways identified with the \code{\link{oranges}} package.
#'
#' @param oranges_res ORANGES data frame
#' @param thr Significance (FDR-corrected q-value) for pathway enrichment.  Defaults to q < 0.05
plot_oranges <- function(oranges_res,thr)
{
  if(missing(thr)) thr <- 0.05
  
  oranges_res %>% 
    dplyr::filter(.,q_val < thr) %>% 
    dplyr::slice(1:25) %>% 
    ggplot(aes(y = -log10(q_val), x = fct_reorder(pathway_name, -log10(q_val)))) + 
    geom_col(alpha = 0.8) + 
    labs(y = "-log10(q value)", x = "Pathway") + 
    coord_flip() + 
    theme_bw() + 
    theme(axis.title = element_text(face = "bold", colour = "black", size = 20),
          axis.text = element_text(colour = "black", size = 12),
          panel.grid = element_blank()) 
  
}

