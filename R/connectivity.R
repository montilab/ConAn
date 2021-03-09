#' @keywords internal
atanh_lower_tri_cor <- function(edat) {
    edat %>%
    stats::cor(method="pearson") %>%
    lower_tri(diag=FALSE) %>%
    remove_na() %>%
    atanh()
}

#' @keywords internal
bg_corrected_atanh_lower_tri_cor <- function(edat, sh) {
    edat %>%
    atanh_lower_tri_cor() %>%
    `*`(sh)
        
}

#' @keywords internal
atanh_lower_tri_erase_mods_cor <- function(edat, mods) {
    edat %>%
    stats::cor(method="pearson") %>%
    erase_mods(mods=mods) %>%
    lower_tri(diag=FALSE) %>%
    remove_na() %>%
    atanh()
}

#' @keywords internal
modular_differential_connectivity <- function(r_edat, t_edat, sh_r, sh_t, mdc_type) {
    
    mc_r <- bg_corrected_atanh_lower_tri_cor(r_edat, sh_r) %>%
            tanh() %>%
            square() %>%
            mean() 

    mc_t <- bg_corrected_atanh_lower_tri_cor(t_edat, sh_t) %>%
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
lower_tri_erase_mods_cor <- function(edat, mods) {
    edat %>%
        stats::cor(method="pearson") %>%
        erase_mods(mods=mods) %>%
        lower_tri(diag=FALSE) %>%
        remove_na()
}

