#' @title Get Connectivity Vectors
#' @description
#' @param mod_genes A character vector of genes within a module
#' @param r_eset Expression set for the reference samples
#' @param t_eset Expression set for the test samples
#' @return A list of connectivity vectors
#'
#' @keywords internal
#'
#' @importFrom magrittr %>%
get_cvs <- function (mod_genes, r_eset, t_eset) {

    # Compute normal adjacency matrix
    r_adj <- r_eset[mod_genes,] %>%
             cor_t_exprs() %>%
             '^'(2)

    # Compute normal adjacency matrix
    t_adj <- t_eset[mod_genes,] %>%
             cor_t_exprs() %>%
             '^'(2)

    # Get upper triangle from adjacency matrix
    cv_r <- get_upper_tri(r_adj, diag=F)
    cv_t <- get_upper_tri(t_adj, diag=F)

    return(cvs = list(cv_r=cv_r, cv_t=cv_t))
}
