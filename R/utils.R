#' @title Get Lower Triangle
#' @description Returns the lower triangle of a matrix in the form of a vector.
#' @param mat A matrix.
#' @param diag Logical. Should the diagonal be included?
#' @return A vector for values.
#'
#' @keywords internal
#'
get_lower_tri <- function(mat, diag) {
  	mat[lower.tri(mat, diag=diag)]
}

#' @keywords internal
remove_na <- function(v) {
    v <- v[!is.na(v)]
    return(v)
}

#' @keywords internal
erase_mods <- function(mat, mod_list) {
    for (i in mod_list) {
      mat[i,i] <- NA
    }
    return(mat)
}