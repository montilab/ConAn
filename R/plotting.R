#' @import magrittr
#' @import ggplot2
#'
#' @export
plot_connectivity <- function(output,N_genes) {
    mapply(function(cv_r, cv_t, mod_name){


        r_name <- output$args$ctrl
        t_name <- output$args$cond

        df <- data.frame(Connectivity = c(cv_r, cv_t),
                         Group = c(rep(r_name, length(cv_r)),
                                   rep(t_name, length(cv_t))))


        if (output$args$mean_correct) {

            cv_r_bg <- unlist(lapply(output$bg$cv_r_bg, mean))
            cv_t_bg <- unlist(lapply(output$bg$cv_t_bg, mean))


            r_bg_name <- paste(r_name, "(BG)")
            t_bg_name <- paste(t_name, "(BG)")

                df <- rbind(df, data.frame(Connectivity = c(cv_r_bg, cv_t_bg),
                                           Group = c(rep(r_bg_name, length(cv_r_bg)),
                                                     rep(t_bg_name, length(cv_t_bg)))))


        }

        # Wes Anderson: Darjeeling
        # Blue Red Teal Orange
        palette <- c("#5BBCD6", "#FF0000", "#00A08A", "#F2AD00")

        df %>%
        ggplot(aes(y=Connectivity, x=Group, fill=Group)) +
        geom_boxplot() +
        theme_minimal() +
        theme(legend.position="none") +
        ylab("Measured Connectivity") +
        ggtitle(mod_name) +
        theme(axis.text.x=element_text(angle=30))+
        scale_fill_manual(values=palette)

    }, output$stat$mods_cv_r, output$stat$mods_cv_t, output$data$mod_names, SIMPLIFY=FALSE)
}

#' @import magrittr
#' @import ggplot2
#'
#' @export
plot_permutations <- function(output) {
    mdc_type <- output$args$mdc_type
    mod_names <- output$data$mod_names
    mdc_permutated <- data.frame(output$iter) %>%
                      set_colnames(mod_names)

    mapply(function(mdc_permuted, mdc_value, mod_name) {

        df <- data.frame(permutations=mdc_permuted)

        if (mdc_type == "fraction") {
            color <- ifelse(mdc_value > 1, "#FF0000", "#5BBCD6")
        }
        if (mdc_type == "difference") {
            color <- ifelse(mdc_value > 0, "#FF0000", "#5BBCD6")
        }

        df %>%
        ggplot(aes(x=permutations)) +
        theme_minimal() +
        geom_density(color="white", fill=color) +
        ggtitle(mod_name) +
        ylab("Probability Density") +
        xlab("Permutated Differential Connectivity") +
        geom_vline(xintercept=mdc_value, linetype="dotted", size=1) +
        theme(axis.title.x=element_text(vjust=-1))

    }, mdc_permutated, output$stat$mods_mdc_adj, mod_names, SIMPLIFY=FALSE)
}
