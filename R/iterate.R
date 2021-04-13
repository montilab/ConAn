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
do_background <- function(iter, c_edat, mods, mean_correct, N_genes=NULL, method) {
	genes <- colnames(c_edat)
	alt_samp <- !is.null(N_genes)

	if(alt_samp) {
		if(N_genes > length(genes)) { stop(paste("N_genes value", N_genes, "is greater than the", length(genes), "number of genes in ExpressionSet object")) }
	}

    if (mean_correct) {
      g_sbst <- if (alt_samp) sample(1:length(genes), N_genes) else 1:length(genes)

      r_m <- c_edat[iter$samples_r, g_sbst]
	  t_m <- c_edat[iter$samples_t, g_sbst]

	  bg_r_cv <- r_m %>%
	      lower_tri_erase_mods_cor(mods=mods, method=method)
	  bg_t_cv <- t_m %>%
	      lower_tri_erase_mods_cor(mods=mods, method=method)

	  iter$bg_r <- bg_r_cv %>%
	      `^`(2) %>%
	      mean
	  iter$bg_t <- bg_t_cv %>%
          `^`(2) %>%
          mean

      # Calculate shrinking factors
	  sh_vec <- get_shrink(bg_r_cv, bg_t_cv, iter$bg_r, iter$bg_t)
	  iter$sh_r <- sh_vec[1]
	  iter$sh_t <- sh_vec[2]

    } else {
        iter$bg_r <- 0
        iter$bg_t <- 0
        iter$sh_r <- 1
        iter$sh_t <- 1
    }

    return(iter)
}

#' @keywords internal
do_differential_connectivity <- function(iter_input, c_edat, mods, mdc_type, method) {

    r_edat <- c_edat[iter_input$samples_r,]
    t_edat <- c_edat[iter_input$samples_t,]

    sh_r <- iter_input$sh_r
    sh_t <- iter_input$sh_t

    mods_mdc <- lapply(mods, function(mod) {
        modular_differential_connectivity(r_edat=r_edat[,mod],
                                          t_edat=t_edat[,mod],
                                          sh_r=sh_r,
                                          sh_t=sh_t,
                                          mdc_type=mdc_type,
                                          method=method)
    })

    iter_output <- list(mods_mdc)
    return(iter_output)
}
