# Get mean connectivity
get_mc <- function (cv) {
  return(mean(cv, na.rm=TRUE))
}
lapply_get_mc <- function (mod_cvs) {
    return(list(mc_r=get_mc(mod_cvs$cv_r),
                mc_t=get_mc(mod_cvs$cv_t)))
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