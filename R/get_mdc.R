#' @title Temporary
#' @description 
#' @param placeholder A placeholder
#' @return A return
#'
#' @keywords internal
#'
#' @importFrom magrittr %>%
#' @importFrom stats ks.test
get_mc <- function (cv) {
  return(mean(cv))
}
lapply_get_mc <- function (mod_cvs) {
    return(list(mc_r=get_mc(mod_cvs$cv_r),
                mc_t=get_mc(mod_cvs$cv_t)))
}

# Get ks score
get_ks <- function (x, y) {
    return(as.numeric(ks.test(x, y)$statistic))
}
lapply_get_ks <- function (mod_cvs) {
    return( get_ks(mod_cvs$cv_r,
                   mod_cvs$cv_t) )
}

# Get module differential connectivity
get_mdc <- function (mc_r, mc_t, bg_r=0, bg_t=0, type="frac") {
  if (type == "frac") {
    mdc <- (mc_t - bg_t) / (mc_r - bg_r)
  }
  if (type == "diff") {
    mdc <- (mc_t - bg_t) - (mc_r - bg_r)
  }
  return(mdc)
}
lapply_get_mdc <- function (mod_mcs, bg_r, bg_t, type) {
  return( get_mdc(mod_mcs$mc_r,
                  mod_mcs$mc_t,
                  bg_r,
                  bg_t,
                  type) )
}

# Get module differential connectivity for each module
get_mods_mdc <- function (mod_list, r_eset, t_eset, mean_correct, mdc_type) {

  if (mean_correct) {
    # Compute adjacency matrix
    r_adj <- cor_t_exprs(r_eset)
    t_adj <- cor_t_exprs(t_eset)

    # For each gene pair in each module
    for (i in mod_list) {
      r_adj[i,i] <- NA
      t_adj[i,i] <- NA
    }

    bg_r <- r_adj %>%
            get_upper_tri(diag=FALSE) %>%
            mean(na.rm=TRUE) %>%
            abs()

    bg_t <- t_adj %>%
            get_upper_tri(diag=FALSE) %>%
            mean(na.rm=TRUE) %>%
            abs()
  } else {
    bg_r <- 0
    bg_t <- 0
  }

  # Output for each module
  mods_cvs <- lapply(mod_list, get_cvs, r_eset, t_eset)
  mods_mcs <- lapply(mods_cvs, lapply_get_mc)
  mods_ks  <- lapply(mods_cvs, lapply_get_ks)
  mods_mdc <- lapply(mods_mcs, lapply_get_mdc, bg_r, bg_t, mdc_type)

  mdc_iter <- list(mods_mdc, mods_ks)
  return(mdc_iter)
}