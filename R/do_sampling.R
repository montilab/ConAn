#' @title Randomly Shuffle Samples
#' @description 
#' @param c_eset Combined expression set control and condition samples
#' @param r_samples Control sample names
#' @param t_samples Condition sample names
#' @param method Sampling method can be either c("bootstrap", "permutation")
#'
#' @return Return shuffled samples for reference and test expression sets
#'
#' @keywords internal
#'
#' @importFrom Biobase sampleNames
do_sampling <- function(c_eset, r_samples, t_samples, method) {
    c_samples <- Biobase::sampleNames(c_eset)
    c_n <- length(c_samples)
    r_n <- length(r_samples)
    t_n <- length(t_samples)

    if (method == "bootstrap") {
        # Sample rerence samples with replacement
        r_samples <- r_samples[sample(1:r_n, r_n, replace=T)]
        # Sample t samples with replacement
        t_samples <- t_samples[sample(1:t_n, t_n, replace=T)]
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