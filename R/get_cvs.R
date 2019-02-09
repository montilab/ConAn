#' @title Get Connectivity Vectors
#' @description 
#' @param mod_genes A character vector of genes within a module
#' @param r_edat A sample x gene expression matrix for the reference samples
#' @param t_edat A sample x gene expression matrix for the test samples
#' @return A list of connectivity vectors
#'
#' @keywords internal
#'
#' @importFrom magrittr %>%
get_cvs <- function (mod_genes, r_edat, t_edat) {

    # Compute normal adjacency matrix
    cv_r <- r_edat[,mod_genes] %>%
            cv()

    # Compute normal adjacency matrix
    cv_t <- t_edat[,mod_genes] %>%
            cv()

    return(cvs = list(cv_r=cv_r, cv_t=cv_t))
}