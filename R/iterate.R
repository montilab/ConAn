#' @title Randomly Shuffle Samples
#' @param c_samples Combined control and condition sample names
#' @param r_samples Control sample names
#' @param t_samples Condition sample names
#' @param method Sampling method can be either c("bootstrap", "permutation")
#'
#' @return Return shuffled samples for reference and test expression sets
#'
#' @keywords internal
do_sampling <- function(iter, c_samples, r_samples, t_samples, method=c("bootstrap", "permutation")) {
    c_n <- length(c_samples)
    r_n <- length(r_samples)
    t_n <- length(t_samples)

    if (method == "bootstrap") {
        r_samples <- r_samples[sample(1:r_n, r_n, replace=TRUE)]
        t_samples <- t_samples[sample(1:t_n, t_n, replace=TRUE)]
    }
    if (method == "permutation") {
        Sample <- sample(1:c_n, c_n)
        r_samples <- c_samples[Sample[1:r_n]]
        t_samples <- c_samples[Sample[(r_n+1):(r_n+t_n)]]
    }

    output <- list()
    output[['samples_r']] <- match(r_samples, c_samples)
    output[['samples_t']] <- match(t_samples, c_samples)

    return(output)
}

#' @keywords internal
do_background <- function(iter, c_edat, mods, mean_correct) {

    if (mean_correct) {
        
        iter[['bg_r']] <- c_edat[iter$samples_r,] %>%
                          atanh_lower_tri_erase_mods_cor(mods=mods) %>%
                          mean()

        iter[['bg_t']] <- c_edat[iter$samples_t,] %>%
                          atanh_lower_tri_erase_mods_cor(mods=mods) %>%
                          mean()
                          
    } else {
        iter[['bg_r']] <- 0
        iter[['bg_t']] <- 0
    }

    return(iter)
}

#' @keywords internal
do_differential_connectivity <- function(iter_input, c_edat, mods, mdc_type) {

    r_edat <- c_edat[iter_input$samples_r,]
    t_edat <- c_edat[iter_input$samples_t,]

    bg_r <- iter_input$bg_r
    bg_t <- iter_input$bg_t

    mods_mdc <- lapply(mods, function(mod) {
        modular_differential_connectivity(r_edat=r_edat[,mod],
                                          t_edat=t_edat[,mod],
                                          bg_r=bg_r,
                                          bg_t=bg_t,
                                          mdc_type=mdc_type)
    })

    iter_output <- list(mods_mdc)
    return(iter_output)
}
