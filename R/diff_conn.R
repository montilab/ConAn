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
                      mean_correct = F,
                      cores = 1,
                      use_gpu = F,
                      mdc_type = c("frac", "diff"),
                      plotting = F,
                      reporting = F,
                      report_dir = "./data") {

    # Grab multi-optional variables
    sim_type <- match.arg(sim_type)
    mdc_type <- match.arg(mdc_type)

    # Create separate expression sets for each condition
    r_eset <- eset[,pData(eset)[,covariate] == ctrl]
    t_eset <- eset[,pData(eset)[,covariate] == cond]

    # Remove genes with no variation within either group
    genes.remove <- genes_novar(r_eset, t_eset)
    r_eset <- r_eset[!(rownames(r_eset) %in% genes.remove),]
    t_eset <- t_eset[!(rownames(t_eset) %in% genes.remove),]
    cat("Removed", length(genes.remove), "genes with no variance in reference/test subsets...\n")

    # Store all data in output variable
    output <- list()
    output$info <- list(mod_names = names(mod_list), # List of module names
                        mod_n = sapply(mod_list, length)) # Size of each module

    # ------------------------------------------
    #    Calculating Background Connectivity
    # ------------------------------------------

    # Calculate background mdc [Expensive]
    # 2x 25000x25000 (e.g. RNA-seq Experiments)
    # ~ 5G each matrix in memory
    # gpuCOR would be good here
    if (mean_correct){
        cat("Calculating background connectivity...\n")

        r_adj <- cor_t_exprs(r_eset)
        t_adj <- cor_t_exprs(t_eset)

        for(i in mod_list) {
            t_adj[i,i] <- NA
            r_adj[i,i] <- NA
        }

        # Background connectivity vector
        cv_r_bg <- r_adj %>%
                   get_upper_tri(diag=F) %>%
                   .[!is.na(.)]

        # Background module connectivity
        mc_r_bg <- cv_r_bg %>%
                   mean() %>%
                   abs()

        # Background connectivity vector
        cv_t_bg <- t_adj %>%
                   get_upper_tri(diag=F) %>%
                   .[!is.na(.)]

        # Background module connectivity
        mc_t_bg <- cv_t_bg %>%
                   mean() %>%
                   abs()

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

    mods_cvs <- lapply(mod_list, get_cvs, r_eset, t_eset)
    mods_mcs <- lapply(mods_cvs, lapply_get_mc)
    mods_mdc_org <- lapply(mods_mcs, lapply_get_mdc, 0, 0, mdc_type)
    mods_mdc_adj <- lapply(mods_mcs, lapply_get_mdc, mc_r_bg, mc_t_bg, mdc_type)
    mods_ks  <- lapply(mods_cvs, lapply_get_ks)

    output$stat <- list(# Reference connvectivity vector for each module
                        mods_cv_r = do.call(rbind, mods_cvs)[,'cv_r'],
                        # Test connvectivity vector for each module
                        mods_cv_t = do.call(rbind, mods_cvs)[,'cv_t'],
                        # Reference module connectivity for each module
                        mods_mc_r = unlist(do.call(rbind, mods_mcs)[,'mc_r']),
                        # Test module connectivity for eacj module
                        mods_mc_t = unlist(do.call(rbind, mods_mcs)[,'mc_t']),
                        # Module differential connectivity for each module
                        mods_mdc_org = unlist(mods_mdc_org),
                        # Adjusted module differential connectivity for each module
                        mods_mdc_adj = unlist(mods_mdc_adj),
                        # Ks score for each module
                        mods_ks = unlist(mods_ks))

    # ------------------------------------------
    #  Permutation-based Significance Testing
    # ------------------------------------------

    cat("Estimating p-values using", iter, "iterations. \n")

    ############################################
    # Start paralellization
    if (!use_gpu) {
        iter_out <- mclapply(seq(iter),
                             do_iter,
                             mod_list = mod_list,
                             r_eset = r_eset,
                             t_eset = t_eset,
                             sim_type = sim_type,
                             mean_correct = mean_correct,
                             mdc_type = mdc_type,
                             mc.cores = cores)
    }
    # End paralellization
    ############################################

    # Results for each iter for each module
    rbind_iter_out <- do.call(rbind, iter_out)

    # For iter for module -> module differential connectivity
    iter_mdc <- rbind_iter_out[,1] %>%
                unlist() %>%
                matrix(nrow=length(names(mod_list))) %>%
                t()

    # For iter for module -> ks statistic
    iter_ks <- rbind_iter_out[,2] %>%
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
        mdc_pval <- pnorm(abs(mdc_stat), mean=0, sd=1, lower.tail=F) * 2
        mdc_fdr <- p.adjust(mdc_pval, method = "BH")

        # Calculate 2 sided p-value for ks
        ks_stdev <- apply(iter_ks, 2, sd)
        ks_stat <- output$stat$mods_ks / ks_stdev
        ks_pval <- pnorm(abs(ks_stat), mean=0, sd=1, lower.tail=F) * 2
        ks_fdr <- p.adjust(ks_pval, method = "BH")

        # Combine stats into a dataframe
        output$iter <- data.frame(mdc_stat = mdc_stat,
                                  mdc_stdev = mdc_stdev,
                                  mdc_pval = mdc_pval,
                                  mdc_fdr = mdc_fdr,
                                  ks_stat = ks_stat,
                                  ks_stdev = ks_stdev,
                                  ks_pval = ks_pval,
                                  ks_fdr = ks_fdr)
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

        # Calculate 2 sided p-value for ks
        iter_ks <- rbind(iter_ks, output$stat$mods_ks) # Prevent p-value = zero
        ks_pval <- apply(iter_ks, 2, function(x){
                            test <- ecdf(x)
                            quant <- test(x[length(x)])
                            return(apply(cbind(quant, 1-quant+(1/length(x))), 1, min)*2)
                        }
        )
        ks_fdr <- p.adjust(ks_pval, method = "BH")

        # Combine stats into a dataframe
        output$iter <- data.frame(mdc_pval = mdc_pval,
                                  mdc_fdr = mdc_fdr,
                                  ks_pval = ks_pval,
                                  ks_fdr = ks_fdr)
    }

    # ------------------------------------------
    #          Plotting and Reporting
    # ------------------------------------------

    if (plotting) {
        cat("Generating plots...\n")
        output$plots <- plot_mdc(output, mean_correct)
    }
    if (reporting) {
        cat("Generating report...\n")
        do_report(output, report_dir)
    }

    cat("Successful finish...\n")
    return(output)
}