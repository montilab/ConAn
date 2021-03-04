#' @keywords internal
shrink_factor <- function(targ, cv) {

    if(targ >= 1 | targ <= 0 ) stop("Target must be between 0 and 1." )
    optimize(opt_shrink, interval = c(0, 1), targ, cv)$minimum

}

#' @keywords internal
opt_shrink <- function(sf, targ, cv) {

    zTrans <- atanh(cv) * sf
    est <- atanh(cv) %>%
        `*`(sf) %>%
        tanh %>%
        `^`(2) %>%
        mean
    abs(targ - est)

}

#' @keywords internal
get_shrink <- function(cv_r, cv_t, mc_r, mc_t) {

    list_cv_bg <- list(cv_r, cv_t)
    mc_vec <- c(mc_r, mc_t)
    mc_min <- min(mc_vec)
    mc_min_log <- mc_vec == mc_min
    sv <- c(1, 1)
    sv[!mc_min_log] <- shrink_factor(mc_vec[mc_min_log], unlist(list_cv_bg[!mc_min_log]))
    return(sv)

}
