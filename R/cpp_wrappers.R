#' @import Rcpp
#' @useDynLib ConAn
#' @export
C_pcor <- function(x) {
    .Call("S_pcor", R_x=x, PACKAGE="ConAn")
}

#' @import Rcpp
#' @useDynLib ConAn
#' @export
C_erase_mods <- function(x, ix) {
    .Call("S_erase_mods", R_x=x,  R_ix=ix, PACKAGE="ConAn")
}

#' @import Rcpp
#' @useDynLib ConAn
#' @export
C_atanh_lower_tri_pcor <- function(x) {
    .Call("S_atanh_lower_tri_pcor", R_x=x, PACKAGE="ConAn")
}

#' @import Rcpp
#' @useDynLib ConAn
#' @export
C_atanh_lower_tri_erase_mods_pcor <- function(x, ix) {
    .Call("S_atanh_lower_tri_erase_mods_pcor", R_x=x, R_ix=ix, PACKAGE="ConAn")
}

#' @import Rcpp
#' @useDynLib ConAn
#' @export
C_mean_atanh_lower_tri_erase_mods_pcor <- function(x, ix) {
    .Call("S_mean_atanh_lower_tri_erase_mods_pcor", R_x=x, R_ix=ix, PACKAGE="ConAn")
}

#' @import Rcpp
#' @useDynLib ConAn
#' @export
C_bg_corrected_atanh_lower_tri_pcor <- function(x, bg) {
    .Call("S_bg_corrected_atanh_lower_tri_pcor", R_x=x, R_bg=bg, PACKAGE="ConAn")
}

#' @import Rcpp
#' @useDynLib ConAn
#' @export
C_mean_bg_corrected_atanh_lower_tri_pcor <- function(x, bg) {
    .Call("S_mean_bg_corrected_atanh_lower_tri_pcor", R_x=x, R_bg=bg, PACKAGE="ConAn")
}

#' @import Rcpp
#' @useDynLib ConAn
#' @export
C_modular_differential_connectivity <- function(xr, xt, bgr, bgt, type) {
    if (type == "frac") {type = 1}
    if (type == "diff") {type = 2}
    .Call("S_modular_differential_connectivity", R_xr=xr, R_xt=xt, R_bgr=bgr, R_bgt=bgt, R_type=type, PACKAGE="ConAn")
}
