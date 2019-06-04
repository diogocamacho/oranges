# ORANGES: Over Representation ANalyses for Gene Expression Signatures

The `ORANGES` package (for **O**ver **R**epresentation **AN**alyses for **G**ene **E**xpression **S**ignatures) is a package designed to do over-representation analyses for an omics signature (transcriptomics or proteomics). This is done using the hypergeometric distribution test (or Fisher's exact test) against a collection of pathway gene sets. `ORANGES` was developed for gene set enrichment of human signatures and the pathway sets were gathered from the [Consensus Pathway DataBase](http://cpdb.molgen.mpg.de). I will describe an example on how to run the ORANGES package using a gene expression signature for non-small cell lung cancer for non-smoking Asian females. This data set can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19804) and can be explored deeply with packages such as `GEOquery` and `limma`.

## Define query set
For the data provided as example here, I have defined 2 groups (tumor and healthy) and determined the differential expression of every gene using the `limma` package. I defined a fold change threshold at log2(fold) > 1 and an FDR corrected p-value threshold p < 0.01 to define the gene expression signature.

```{r define_signature}
query <- lung_res %>%
  tibble::add_column(., entrez_id = genes_lung$ENTREZID, .before = 1) %>%
  tibble::add_column(., symbol = genes_lung$ENTREZID, .before = 2) %>%
  dplyr::filter(., abs(logFC) > 1, adj.P.Val < 0.01)
```

Below I illustrate what this looks like in a volcano plot, where the genes in red are those that are inlcuded in the signature.

```{r volcano}
lung_res %>%
  ggplot() +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val)), alpha = 0.5) +
  geom_point(data = subset(lung_res, abs(logFC) > 1 & adj.P.Val < 0.01), aes(x = logFC, y = -log10(adj.P.Val)), alpha = 0.5, color = "red") +
  labs(x = "log2(fold change)", y = "-log10(adjusted p-value)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))
```


## Perform enrichment
With the signature defined, we will work our way through the `ORANGES` functions.

### Subset pathways
The first step in `ORANGES` is to define which of the genes in our data are present in the pathway sets and, of these, which ones are genes from our signature in those pathways as well. We do this with the `subset_pathways` and `extract_pathways` functions:

```{r subset_pathways}
wp <- subset_pathways(query_entrez = query$entrez_id, universe_entrez = genes_lung$ENTREZID)
eps <- extract_pathways(wp)
```

### Genes per pathway
Those two steps above allowed us to have a count for how many genes in our data are present in each pathway and how many from our signature are present in those pathways. We can use the `genes_pathway` function to extract the gene symbols for these:

```{r genes_pathway}
gid <- genes_pathway(genes_query = eps$genes_query, pathways_query = eps$pathways_query)
```

Visualizing them we get:

```{r gene_symbols_pathways}
head(gid)
```


### Calculate enrichment
We can now calculate the pathway enrichment using Fisher's exact test:

```{r calculate_enrichment}
enr <- calculate_enrichment(genes_query = eps$genes_query,
                              pathways_universe = eps$pathways_universe,
                              pathways_query = eps$pathways_query)
```

which provides us with a data frame where each row corresponds to a given pathway for a given data source, and where the nominal and corrected p-values for the enrichment are computed.

```{r enr_table}
head(enr)
```

This is hard to read and so next we will generate a results data frame.

### Results data frame
Taking the output of the functions above, we can construct a results data frame with the `enrichment_results` function:

```{r}
res <- enrichment_results(genes_pathway = gid,
                            number_genes = enr$number_genes,
                            gene_ratio = enr$gene_ratio,
                            nominal_p = enr$nominal_p,
                            p_adj = enr$padj)
```

which is dramatically easier to understand than the previous table, as illustrated here:

```{r res_table}
head(res)
```

We can finally produce a plot of the enrichment using `ggplot2`, as you'll see next.

### Plotting
Plotting of the enrichment results can be easily done with the `plot_oranges` function:

```{r plotting}
plot_oranges(res)
```
