#' @title No-Variance Genes
#' @description Finds genes without variance across either control or condition samples
#' @param r_eset Expression set of control samples
#' @param t_eset Expression set of condition samples
#' @return A vector of gene names with no variance
#'
#' @keywords internal
#'
#' @importFrom Biobase exprs
genes_novar <- function(r_eset, t_eset) {
    genes <- rownames(r_eset)
    genes.remove <- c()
    genes.remove <- c(genes.remove, genes[apply(exprs(r_eset), 1, var) == 0])
    genes.remove <- c(genes.remove, genes[apply(exprs(t_eset), 1, var) == 0])
    return(genes.remove)
}