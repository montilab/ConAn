#' @title Get Lower Triangle
#' @description Returns the lower triangle of a matrix in the form of a vector.
#' @param mat A matrix.
#' @param diag Logical. Should the diagonal be included?
#' @return A vector for values.
#'
#' @keywords internal
lower_tri <- function(mat, diag) {
  	mat[lower.tri(mat, diag=diag)]
}

#' @keywords internal
remove_na <- function(v) {
    v[!is.na(v)]
}

#' @keywords internal
subtract_bg <- function(mat, bg) {
    mat-bg
}

#' @keywords internal
square <- function(x) {
    x^2
}

#' @keywords internal
erase_mods <- function(mat, mods) {
    for (i in mods) {
      mat[i,i] <- NA
    }
    return(mat)
}