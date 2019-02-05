#' @title Randomly Shuffle Samples
#' @description 
#' @param c_samples Combined control and condition sample names
#' @param r_samples Control sample names
#' @param t_samples Condition sample names
#' @param method Sampling method can be either c("bootstrap", "permutation")
#'
#' @return Return shuffled samples for reference and test expression sets
#'
#' @keywords internal
do_sampling <- function(c_samples, r_samples, t_samples, method) {
    c_n <- length(c_samples)
    r_n <- length(r_samples)
    t_n <- length(t_samples)

    if (method == "bootstrap") {
        # Sample rerence samples with replacement
        r_samples <- r_samples[sample(1:r_n, r_n, replace=TRUE)]
        # Sample t samples with replacement
        t_samples <- t_samples[sample(1:t_n, t_n, replace=TRUE)]
    }
    if (method == "permutation") {
        #Sample rerence samples with replacement
        Sample <- sample(1:c_n, c_n)
        r_samples <- c_samples[Sample[1:r_n]]
        t_samples <- c_samples[Sample[(r_n+1):(r_n+t_n)]]
    }

    sample_order <- list()
    sample_order[['r']] <- match(r_samples, c_samples)
    sample_order[['t']] <- match(t_samples, c_samples)

    return(sample_order)
}