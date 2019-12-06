#' @title Randomly Shuffle Samples
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

    bg_r <- c_edat[iter_input$samples_r,] %>%
            stats::cor() %>%
            erase_mod_list(mod_list=mod_list) %>%
            lower_tri(diag=FALSE) %>%
            remove_na() %>%
            atanh() %>%
            mean()

    bg_t <- c_edat[iter_input$samples_t,] %>%
            stats::cor() %>%
            erase_mod_list(mod_list=mod_list) %>%
            lower_tri(diag=FALSE) %>%
            remove_na() %>%
            atanh() %>%
            mean()

    iter_output <- list()
    iter_output[['samples_r']] <- iter_input$samples_r
    iter_output[['samples_t']] <- iter_input$samples_t
    iter_output[['bg_r']] <- bg_r
    iter_output[['bg_t']] <- bg_t

    return(iter_output)
}

skip_background <- function(iter_input) {

    iter_output <- list()
    iter_output[['samples_r']] <- iter_input$samples_r
    iter_output[['samples_t']] <- iter_input$samples_t
    iter_output[['bg_r']] <- 0
    iter_output[['bg_t']] <- 0

    return(iter_output)
}

# Get module differential connectivity for each module
do_differential_connectivity <- function(iter_input, c_edat, mod_list, mdc_type) {

    r_edat <- c_edat[iter_input$samples_r,]
    t_edat <- c_edat[iter_input$samples_t,]

    bg_r <- iter_input$bg_r
    bg_t <- iter_input$bg_t

    mods_mdc <- lapply(mod_list, function(mod_genes) {

        bccv.r <- r_edat[,mod_genes] %>%
                  stats::cor() %>%
                  lower_tri(diag=FALSE) %>%
                  remove_na() %>%
                  atanh() %>%
                  subtract_bg(bg_r)    
            
        bccv.t <- t_edat[,mod_genes] %>%
                  stats::cor() %>%
                  lower_tri(diag=FALSE) %>%
                  remove_na() %>%
                  atanh() %>%
                  subtract_bg(bg_t)  

        mc.r <- mean( tanh(bccv.r)^2 )
        mc.t <- mean( tanh(bccv.t)^2 )

        if (mdc_type == "frac") { 
            return(mc.t / mc.r)
        }
        if (mdc_type == "diff") {
            return(mc.t - mc.r)
        }
    })

    iter_output <- list(mods_mdc)
    return(iter_output)
}
