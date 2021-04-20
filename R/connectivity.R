#' @keywords internal
atanh_lower_tri_cor <- function(edat, corr_func) {
    print("atanh_lower_tri_cor")
    edat %>%
    corr_func() %>%
    lower_tri(diag=FALSE) %>%
    remove_na() %>%
    atanh()
}

#' @keywords internal
bg_corrected_atanh_lower_tri_cor <- function(edat, sh, corr_func) {
    print("bg_corrected_atanh_lower_tri_cor")
    edat %>%
    atanh_lower_tri_cor(corr_func=corr_func) %>%
    `*`(sh)
        
}

#' @keywords internal
atanh_lower_tri_erase_mods_cor <- function(edat, mods, corr_func) {
    print("atanh_lower_tri_erase_mods_cor")
    edat %>%
    corr_func() %>%
    erase_mods(mods=mods) %>%
    lower_tri(diag=FALSE) %>%
    remove_na() %>%
    atanh()
}

#' @keywords internal
modular_differential_connectivity <- function(r_edat, t_edat, sh_r, sh_t, mdc_type, corr_func) {
    print("modular_differential_connectivity")
    
    mc_r <- bg_corrected_atanh_lower_tri_cor(r_edat, sh_r, corr_func) %>%
            tanh() %>%
            square() %>%
            mean() 

    mc_t <- bg_corrected_atanh_lower_tri_cor(t_edat, sh_t, corr_func) %>%
            tanh() %>%
            square() %>%
            mean()

    if (mdc_type == "fraction") { 
        return(mc_t / mc_r)
    }
    if (mdc_type == "difference") {
        return(mc_t - mc_r)
    }
}

#' @keywords internal
lower_tri_erase_mods_cor <- function(edat, mods, corr_func) {
    print("lower_tri_erase_mods_cor")
    edat %>%
        corr_func() %>%
        erase_mods(mods=mods) %>%
        lower_tri(diag=FALSE) %>%
        remove_na()
}

