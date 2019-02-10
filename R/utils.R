#' @title Get Lower Triangle
#' @description Returns the lower triangle of a matrix in the form of a vector.
#' @param mat A matrix.
#' @param diag Logical. Should the diagonal be included?
#' @return A vector for values.
#'
#' @keywords internal
get_lower_tri <- function(mat, diag) {
  	mat[lower.tri(mat, diag=diag)]
}

#' @keywords internal
remove_na <- function(v) {
    v[!is.na(v)]
}

#' @keywords internal
erase_mods <- function(mat, mod_list) {
    for (i in mod_list) {
      mat[i,i] <- NA
    }
    return(mat)
}

#' @importFrom dplyr bind_rows
#' @keywords internal
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