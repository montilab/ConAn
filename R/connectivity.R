modlist.to.matzindex <- function(modlist, genes) {
    modlist.index <- lapply(modlist, function(x) {
        sapply(x, function(y) {
            match(y, genes)
        })
    })
 
    # Create pairwise combinations for each module
    # Bind into a single dataframe
    index.pairs <- modlist.index %>%
                   lapply(function(x) expand.grid(x, x)) %>%
                   dplyr::bind_rows()
            
    # Matrix-like dimension names
    colnames(index.pairs) <- c("i", "j")
    
    # Zero-indexed
    index.pairs <- index.pairs-1
    
    # Matrix is represented as a zero-indexed vector 
    mat.size <- length(genes)
    mat.zindex <- index.pairs$i * mat.size + index.pairs$j

    return(mat.zindex)  
}

#' @import Rcpp
#' @useDynLib ConAn
pcor <- function(x) {
    .Call("S_pcor", R_x=x, PACKAGE="ConAn")
}

#' @import Rcpp
#' @useDynLib ConAn
cv <- function(x) {
    .Call("S_cv", R_x=x, PACKAGE="ConAn")
}

#' @import Rcpp
#' @useDynLib ConAn
bgcv <- function(x, ix) {
    v <- .Call("S_bgcv", R_x=x, R_ix=ix, PACKAGE="ConAn")
    v <- v[v != 0]
    return(v)
}

#' @import Rcpp
#' @useDynLib ConAn
mc <- function(x) {
    .Call("S_mc", R_x=x, PACKAGE="ConAn")
}

#' @import Rcpp
#' @useDynLib ConAn
bgmc <- function(x, ix) {
    .Call("S_bgmc", R_x=x, R_ix=ix, PACKAGE="ConAn")

}