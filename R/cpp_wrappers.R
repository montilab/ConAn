#' @import Rcpp
#' @useDynLib ConAn
#' @export
pcor <- function(x) {
    .Call("S_pcor", R_x=x, PACKAGE="ConAn")
}

#' @import Rcpp
#' @useDynLib ConAn
#' @export
erase_vals <- function(x, ix) {
    .Call("S_erase_vals", R_x=x,  R_ix=ix, PACKAGE="ConAn")
}

#' @import Rcpp
#' @useDynLib ConAn
#' @export
cv <- function(x) {
    .Call("S_cv", R_x=x, PACKAGE="ConAn")
}

#' @import Rcpp
#' @useDynLib ConAn
#' @export
bgcv <- function(x, ix) {
    v <- .Call("S_bgcv", R_x=x, R_ix=ix, PACKAGE="ConAn")
    v <- v[v != 0]
    return(v)
}

#' @import Rcpp
#' @useDynLib ConAn
#' @export
mc <- function(x) {
    .Call("S_mc", R_x=x, PACKAGE="ConAn")
}

#' @import Rcpp
#' @useDynLib ConAn
#' @export
bgmc <- function(x, ix) {
    .Call("S_bgmc", R_x=x, R_ix=ix, PACKAGE="ConAn")

}