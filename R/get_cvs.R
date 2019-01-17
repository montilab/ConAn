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
             cor_t_exprs()

    # Compute normal adjacency matrix
    t_adj <- t_eset[mod_genes,] %>%
             cor_t_exprs()

    # Get upper triangle from adjacency matrix
    cv_r <- r_adj %>%
            get_upper_tri(diag=FALSE) %>%
            atanh()

    cv_t <- t_adj %>%
            get_upper_tri(diag=FALSE) %>%
            atanh()

    return(cvs = list(cv_r=cv_r, cv_t=cv_t))
}
