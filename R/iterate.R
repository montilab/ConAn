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
do_sampling <- function(iter_input, c_samples, r_samples, t_samples, method=c("bootstrap", "permutation")) {
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

    iter_output <- list()
    iter_output[['samples_r']] <- match(r_samples, c_samples)
    iter_output[['samples_t']] <- match(t_samples, c_samples)

    return(iter_output)
}

do_background <- function(iter_input, c_edat, mod_list) {

    # Compute reference background
    bg_r <- c_edat[iter_input$samples_r,] %>%
            cor() %>%
            erase_mods(mod_list=mod_list) %>%
            get_upper_tri(diag=FALSE) %>%
            abs() %>%
            mean(na.rm=TRUE)

    # Compute test background
    bg_t <- c_edat[iter_input$samples_t,] %>%
            cor() %>%
            erase_mods(mod_list=mod_list) %>%
            get_upper_tri(diag=FALSE) %>%
            abs() %>%
            mean(na.rm=TRUE)

    iter_output <- list()
    iter_output[['samples_r']] <- iter_input$samples_r
    iter_output[['samples_t']] <- iter_input$samples_t
    iter_output[['bg_r']] <- bg_r
    iter_output[['bg_t']] <- bg_t

    return(iter_output)
}

skip_background <- function(iter_input) {

    iter_output <- list()
    iter_output[['samples_r']] <- iter_input$r
    iter_output[['samples_t']] <- iter_input$t
    iter_output[['bg_r']] <- 0
    iter_output[['bg_t']] <- 0

    return(iter_output)
}

# Get module differential connectivity for each module
get_mods_mdc <- function (iter_input, c_edat, mod_list, mdc_type) {

    r_edat <- c_edat[iter_input$samples_r,]
    t_edat <- c_edat[iter_input$samples_t,]

    bg_r <- iter_input$bg_r
    bg_t <- iter_input$bg_t

    # Output for each module
    mods_cvs <- lapply(mod_list, get_cvs, r_edat, t_edat)
    mods_mcs <- lapply(mods_cvs, lapply_get_mc)
    mods_ks  <- lapply(mods_cvs, lapply_get_ks)
    mods_mdc <- lapply(mods_mcs, lapply_get_mdc, bg_r, bg_t, mdc_type)

    iter_output <- list(mods_mdc, mods_ks)
    return(iter_output)
}
