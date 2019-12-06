#' @title Plot Module Differential Connectivity
#' @param output Output of diff_conn
#' @param mean_correct Logical. Was mean correct used in diff_conn?
#' @return A list of ggplots
#'
#' @import ggplot2
#'
#' @export
do_plot <- function(output, mean_correct) {
    # Extract module connectivity vectors
    cv_r <- output$stat$mods_cv_r
    cv_t <- output$stat$mods_cv_t
    conn_list <- list()
    for (i in 1:length(cv_r)) {
        conn_list[[i]] <- list(cv_r[[i]], cv_t[[i]])
        names(conn_list[[i]]) <- rep(names(cv_r)[i], 2)
    }

    # Plot PDF
    plot_pdf <- lapply(conn_list, function(x){

    cv_r <- unlist(x[[1]])
    cv_t <- unlist(x[[2]])

    df <- data.frame(Connectivity = c(cv_r, cv_t),
                     Group = c(rep("Control", length(cv_r)),
                               rep("Test", length(cv_t))))


    df$Group <- factor(df$Group, levels=c("Control", "Test", "Control Background", "Test Background"))

    if (mean_correct) {

        cv_r_bg <- output$bg$cv_r_bg
        cv_t_bg <- output$bg$cv_t_bg

        df <- rbind(df, data.frame(Connectivity = c(cv_r_bg, cv_t_bg),
                                 Group = c(rep("Control Background", length(cv_r_bg)),
                                           rep("Test Background", length(cv_t_bg)))))
    }

    p <- ggplot(df, aes(x = Connectivity, fill=Group)) +
                stat_density(alpha=0.5, position="identity") +
                scale_y_continuous(name = "Density") +
                ggtitle(names(x)[1])

    return(p)
    })

    # Plot CDF
    plot_cdf <- lapply(conn_list, function(x){

    cv_r <- unlist(x[[1]])
    cv_t <- unlist(x[[2]])

    df <- data.frame(Connectivity <- c(cv_r, cv_t),
                     Group = c(rep("Control", length(cv_r)),
                               rep("Test", length(cv_r))))

    p <- ggplot(df, aes(x = Connectivity, colour=Group)) +
                stat_ecdf(geom = "step") +
                scale_y_continuous(name = "CDF") +
                ggtitle(names(x)[1])
    return(p)
    })

    names(plot_pdf) <- names(plot_cdf) <- output$info$mod_names

    return(list(plot_pdf = plot_pdf,
                plot_cdf = plot_cdf))
}