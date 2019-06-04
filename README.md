# ORANGES: Over Representation ANalyses for Gene Expression Signatures

The `ORANGES` package is a package designed to do over-representation analyses for an omics signature (transcriptomics or proteomics). This is done using the hypergeometric distribution test (or Fisher's exact test) against a collection of pathway gene sets. `ORANGES` was developed for gene set enrichment of human signatures and the pathway sets were gathered from the [Consensus Pathway DataBase](http://cpdb.molgen.mpg.de).

## Installing
To install the package do, in R:

```{r}
devtools::install_github("diogocamacho/oranges")
```

## Running ORANGES
The easiest way to run the algorithm is to use the `oranges` wrapper. In R, do:

```{r}
enrichment <- oranges(query_entrez, universe_entrez)
```

where `query_entrez` are the EntrezIDs for the query of interest, and `universe_entrez` are the EntrezIDs of the background (e.g., the total set of genes identified in the platform of interest). The output of the wrapper is a tidy data frame where columns are source of the pathway set, name of pathway, gene ratio, and nominal and corrected p-values for the test statistic.

## Vignette
Please see the package vignette for a detailed example:

```{r}
browseVignette("oranges")
```

## Contact
For suggestions, comments, additions, please do not hesitate to [email me](mailto:diogo.camacho@wyss.harvard.edu)
