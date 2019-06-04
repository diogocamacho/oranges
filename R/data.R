#' Pathway sets for ORANGES
#' 
#' Human pathway sets collected from the Consensus Pathway Database
#' 
#' @format A matrix, a data frame, and a character vector:
#' \describe{
#'   \item A matrix (P) with the element assignments for a given gene in a given pathway
#'   \item A data frame (pathway_info) with 3 columns:
#'     \item Pathway source (KEGG, Reactome, etc)
#'     \item Pathway name
#'     \item An internal id reference for the pathway sources 
#'   \item A vector of gene symbols for all gene EntrezIDs in the matrix P.
#' }
#' 
"pathway_sets"
