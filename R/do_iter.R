#' @title Iteration of Random Sampling
#' @description Performs an iteration random sampling and calculating differential connectivity
#' @param n An unused placeholder variable
#' @param mod_list A list of modules of gene names
#' @param r_eset Expression set of control samples
#' @param t_eset Expression set of condition samples
#' @param sim_type Simulation type can be either c("bootstrap", "permutation")
#' @param mean_correct Correct with background connectivity
#' @param mdc_type Method for calculating difference in connectivity can be either c("frac", "diff")
#' 
#' @return Connectivity statistics for the iteration
#'
#' @keywords internal
#'
#' @importFrom BiocGenerics combine
#' @importFrom Biobase sampleNames
do_iter <- function(n, mod_list, r_eset, t_eset, sim_type, mean_correct, mdc_type) {
    c_eset <- BiocGenerics::combine(r_eset, t_eset)
    r_samples <- Biobase::sampleNames(r_eset)
    t_samples <- Biobase::sampleNames(t_eset)
    sample_order <- do_sampling(c_eset, r_samples, t_samples, sim_type)
    r_eset <- c_eset[,sample_order$r]
    t_eset <- c_eset[,sample_order$t]
    iter_out <- get_mods_mdc(mod_list, r_eset, t_eset, mean_correct, mdc_type)
}