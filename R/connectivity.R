#' @import bigcor
#' @import ffbase

#' @keywords internal
atanh_lower_tri_cor <- function(edat, bigcor_on) {
    if (bigcor_on){
        edat %>%
            bigcor(size = dim(edat)[2]) %>%
            lower_tri(diag=FALSE,bigcor_on) %>%
            remove_na() %>%
            atanh()
    } else {
        edat %>%
            stats::cor(method="pearson") %>%
            lower_tri(diag=FALSE,bigcor_on) %>%
            remove_na() %>%
            atanh()
    }
}

#' @keywords internal
bg_corrected_atanh_lower_tri_cor <- function(edat, bg, bigcor_on) {
    edat %>%
    atanh_lower_tri_cor(bigcor_on) %>%
    subtract_bg(bg)
}

#' @keywords internal
atanh_lower_tri_erase_mods_cor <- function(edat, mods, bigcor_on) {
    if (bigcor_on){
        edat %>%
            bigcor(size = dim(edat)[2]) %>%
            erase_mods(mods=mods) %>%
            lower_tri(diag=FALSE, bigcor_on) %>%
            remove_na() %>%
            atanh()
    } else {
        edat %>%
            stats::cor(method="pearson") %>%
            erase_mods(mods=mods) %>%
            lower_tri(diag=FALSE,bigcor_on) %>%
            remove_na() %>%
            atanh()
    }
    
}

#' @keywords internal
modular_differential_connectivity <- function(r_edat, t_edat, bg_r, bg_t, mdc_type) {
    
    mc_r <- bg_corrected_atanh_lower_tri_cor(r_edat, bg_r) %>%
            tanh() %>%
            square() %>%
            mean(na.rm=TRUE) 

    mc_t <- bg_corrected_atanh_lower_tri_cor(t_edat, bg_t) %>%
            tanh() %>%
            square() %>%
            mean(na.rm=TRUE)

    if (mdc_type == "fraction") { 
        return(mc_t / mc_r)
    }
    if (mdc_type == "difference") {
        return(mc_t - mc_r)
    }
}
