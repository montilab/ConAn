#' @title Differential Connectivity
#' @description Calculates the differential connectivity between two groups of samples across one or more modules of genes.
#' @param eset An expression set object
#' @param mod_list A list of character vectors
#' @param covariate A string corresponding to column in pData distinguishing sample groups
#' @param ctrl A string defining control samples
#' @param cond A string defining condition samples
#' @param sim_type Simulation type can be either c("bootstrap", "permutation")
#' @param iter Number of sampling iterations during significance testing
#' @param mean_correct Correct with background connectivity
#' @param cores Number of CPU cores available for parallelization
#' @param use_gpu Use of GPU computing
#' @param mdc_type Method for calculating difference in connectivity can be either c("frac", "diff")
#' @param reporting Generate a markdown report for analysis
#' @param report_dir Directory where report is generated
#' 
#' @return A list of statistics and plots resulting from the analysis
#'
#' @importFrom Biobase pData
#' @importFrom magrittr %>%
#' @importFrom stats pnorm
#' @importFrom stats sd
#' @importFrom stats ecdf
#' @importFrom stats p.adjust
#' @importFrom parallel mclapply
#'
#' @export
diff_conn <- function(eset,
                      mod_list,
                      covariate,
                      ctrl,
                      cond,
                      sim_type = c("bootstrap", "permutation"),
                      iter = 5,
                      mean_correct = FALSE,
                      cores = 1,
                      use_gpu = FALSE,
                      mdc_type = c("frac", "diff"),
                      plotting = FALSE,
                      reporting = FALSE,
                      report_dir = "./data") {

    # Grab multi-optional variables
    sim_type <- match.arg(sim_type)
    mdc_type <- match.arg(mdc_type)

    pdat <- pData(eset)  
    c_samples <- rownames(pdat)
    r_samples <- c_samples[(pdat[,covariate] == ctrl)]
    t_samples <- c_samples[(pdat[,covariate] == cond)]

    # Extract and format expression matrix
    c_edat <- t(Biobase::exprs(eset)) # A sample x gene expression matrix
    r_edat <- c_edat[r_samples,]
    t_edat <- c_edat[t_samples,]

    # Genes
    genes <- colnames(c_edat)

    # Store all data in output variable
    output <- list()
    output$info <- list(mod_names = names(mod_list), # List of module names
                        mod_n = sapply(mod_list, length)) # Size of each module

    # ------------------------------------------
    #    Calculating Background Connectivity
    # ------------------------------------------

    # Calculate background mdc
    if (mean_correct){
        cat("Calculating background connectivity...\n")

        mat.zindex <- modlist.to.matzindex(mod_list, genes)

        # Background connectivity vector
        cv_r_bg <- atanh_lower_tri_erase_mods_pcor(r_edat, mat.zindex)

        # Background module connectivity
        mc_r_bg <- mean(cv_r_bg)

        # Background connectivity vector
        cv_t_bg <- atanh_lower_tri_erase_mods_pcor(t_edat, mat.zindex)

        # Background module connectivity
        mc_t_bg <- mean(cv_t_bg)

        # Storage for background statistics
        output$bg <- list(cv_r_bg=cv_r_bg,
                          mc_r_bg=mc_r_bg,
                          cv_t_bg=cv_t_bg,
                          mc_t_bg=mc_t_bg)
    } else {
        mc_r_bg = 0
        mc_t_bg = 0
    }

    # ------------------------------------------
    #    Calculating Module Connectivity
    # ------------------------------------------

    cat("Calculating module differential connectivity for", length(names(mod_list)), "clusters...\n")

    # Lambda helper functions
    l_cvs <- function (mod_genes, r_edat, t_edat) {
        cv_r <- r_edat[,mod_genes] %>% bg_corrected_atanh_lower_tri_pcor(mc_r_bg)
        cv_t <- t_edat[,mod_genes] %>% bg_corrected_atanh_lower_tri_pcor(mc_t_bg)
        return(cvs = list(cv_r=cv_r, cv_t=cv_t))
    }

    # Background-corrected modular connectivity vectors
    mods_cvs <- lapply(mod_list, l_cvs, r_edat, t_edat)

    # Background-corrected modular connectivity
    mods_mcs <- lapply(mods_cvs, function(x) {
      return(list(mc_r = mean( tanh( x$cv_r )^2 ),
                  mc_t = mean( tanh( x$cv_t )^2 )))
    })
    lapply_get_mdc <- function (mods_mcs) {
        if (mdc_type == "frac") { return(mods_mcs$mc_t / mods_mcs$mc_r) }
        if (mdc_type == "diff") { return(mods_mcs$mc_t - mods_mcs$mc_r) }
    }
    # Adjusted module differential connectivity 
    mods_mdc_adj <- lapply(mods_mcs, lapply_get_mdc)

    output$stat <- list(# Reference connvectivity vector for each module
                        mods_cv_r = do.call(rbind, mods_cvs)[,'cv_r'],
                        # Test connvectivity vector for each module
                        mods_cv_t = do.call(rbind, mods_cvs)[,'cv_t'],
                        # Reference module connectivity for each module
                        mods_mc_r = unlist(do.call(rbind, mods_mcs)[,'mc_r']),
                        # Test module connectivity for eacj module
                        mods_mc_t = unlist(do.call(rbind, mods_mcs)[,'mc_t']),
                        # Adjusted module differential connectivity for each module
                        mods_mdc_adj = unlist(mods_mdc_adj))

    # ------------------------------------------
    #  Permutation-based Significance Testing
    # ------------------------------------------

    cat("Estimating p-values using", iter, "iterations. \n")
    
    ############################################
    # Start paralellization
    #
    # 1.
    # A list of randomly shuffled groups of samples
    iter_sampling <- mclapply(seq(iter),
                              do_sampling,
                              c_samples = c_samples,
                              r_samples = r_samples,
                              t_samples = t_samples,
                              method = sim_type,
                              mc.cores = cores)

    # 2.
    # Background connectivity for each iteration
    if (mean_correct) {
        iter_background <- mclapply(iter_sampling, do_background, c_edat = c_edat,  mat.zindex = mat.zindex, mc.cores = cores)
    } else {
        iter_background <- mclapply(iter_sampling, skip_background, mc.cores = cores) 
    }
    
    # 3.
    # Calculate differential module connectivity
    iter_out <- mclapply(iter_background,
                       do_differential_connectivity,
                       c_edat = c_edat,
                       mod_list = mod_list,
                       mdc_type = mdc_type,
                       mc.cores = cores)
    #
    #
    # End paralellization
    ############################################

    # Results for each iter for each module
    rbind_iter_out <- do.call(rbind, iter_out)

    # For iter for module -> module differential connectivity
    iter_mdc <- rbind_iter_out[,1] %>%
                unlist() %>%
                matrix(nrow=length(names(mod_list))) %>%
                t()

    # ------------------------------------------
    #           Calculate Significance
    # ------------------------------------------

    if (sim_type == "bootstrap") {

        # Calculate 2 sided p-value for mdc
        mdc_stdev <- apply(iter_mdc, 2, sd)
        m1 <- ifelse(mdc_type == "frac", 1, 0)
        mdc_stat <- (output$stat$mods_mdc_adj - m1) / mdc_stdev
        mdc_pval <- pnorm(abs(mdc_stat), mean=0, sd=1, lower.tail=FALSE) * 2
        mdc_fdr <- p.adjust(mdc_pval, method = "BH")

        # Combine stats into a dataframe
        output$iter <- data.frame(mdc_stat = mdc_stat,
                                  mdc_stdev = mdc_stdev,
                                  mdc_pval = mdc_pval,
                                  mdc_fdr = mdc_fdr)
    }

    if (sim_type == "permutation") {

        # Calculate 2 sided p-value for mdc
        iter_mdc <- rbind(iter_mdc, output$stat$mods_mdc_adj) # Prevent p-value = zero
        mdc_pval <- apply(iter_mdc, 2, function (x) {
                             test <- ecdf(x)
                             quant <- test(x[length(x)])
                             return(apply(cbind(quant, 1-quant+(1/length(x))), 1, min)*2)
                          }
                    )
        mdc_fdr <- p.adjust(mdc_pval, method = "BH")

        # Combine stats into a dataframe
        output$iter <- data.frame(mdc_pval = mdc_pval,
                                  mdc_fdr = mdc_fdr)
    }

    # ------------------------------------------
    #          Plotting and Reporting
    # ------------------------------------------

    if (plotting) {
        cat("Generating plots...\n")
        output$plots <- do_plot(output, mean_correct)
    }
    if (reporting) {
        cat("Generating report...\n")
        do_report(output, report_dir)
    }

    cat("Successful finish...\n")
    return(output)
}