#' @title Expression Correlation
#' @description Calculates the correlation of genes across samples from an expression set object.
#' @param eset An expression set object
#' @return A single numeric value
#'
#' @keywords internal
#'
#' @importFrom Biobase exprs
#' @importFrom magrittr %>%
cor_t_exprs <- function(eset) {
    eset %>%
    Biobase::exprs() %>%
    t() %>%
    cor()
}