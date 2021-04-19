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
multiply <- function(x, m) {
    x*m
}

#' @keywords internal
erase_mods <- function(mat, mods) {
    genes <- rownames(mat)
    for (i in mods) {
      s <- i[i %in% genes]
      mat[s,s] <- NA
    }
    return(mat)
}

#' Format a string using placeholders
#'
#' @param string A an unformatted string with placeholders
#' @param ... Variables to format placeholders with
#' @return A formatted string
#' 
#' @examples
#' \dontrun{
#' format_str("Format with {1} and {2}", "x", "y")
#' }
#'
#' @keywords internal
format_str <- function(string, ...) {
    args <- list(...)
    for (i in 1:length(args)) {
        pattern <- paste("\\{", i, "}", sep="")
        replacement <- args[[i]]
        string <- gsub(pattern, replacement, string)
    }
    return(string)
}
