#' @title Get Upper Triangle
#' @description Returns the upper triangle of a matrix in the form of a vector.
#' @param mat A matrix.
#' @param diag Logical. Should the diagonal be included?
#' @return A vector for values.
#'
#' @keywords internal
#'
get_upper_tri <- function(mat, diag) {
  	mat[upper.tri(mat, diag=diag)]
}

#' @keywords internal
#'
erase_mods <- function(mat, mod_list) {
    for (i in mod_list) {
      mat[i,i] <- NA
    }
    return(mat)
}