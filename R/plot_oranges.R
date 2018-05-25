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

