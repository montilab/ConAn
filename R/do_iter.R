#' @title Iteration of Random Sampling
#' @description Performs an iteration random sampling and calculating differential connectivity
#' @param n An unused placeholder variable
#' @param mod_list A list of modules of gene names
#' @param c_edat A sample x gene expression matrix
#' @param r_samples Reference samples
#' @param t_samples Test samples
#' @param sim_type Simulation type can be either c("bootstrap", "permutation")
#' @param mean_correct Correct with background connectivity
#' @param mdc_type Method for calculating difference in connectivity can be either c("frac", "diff")
#' 
#' @return Connectivity statistics for the iteration
#'
#' @keywords internal
do_iter <- function(n, mod_list, c_edat, r_samples, t_samples, sim_type, mean_correct, mdc_type) {
    c_samples <- rownames(c_edat)
    
    sample_order <- do_sampling(c_samples, r_samples, t_samples, sim_type)

    r_edat <- c_edat[sample_order$r,]
    t_edat <- c_edat[sample_order$t,]

    iter_out <- get_mods_mdc(mod_list, r_edat, t_edat, mean_correct, mdc_type)
}