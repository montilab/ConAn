#' @import magrittr
#' @import ggplot2
#'
#' @export
plot_connectivity <- function(output) {
 mapply(function(cv_r, cv_t, mod_name){
    
    r_name <- output$args$ctrl
    t_name <- output$args$cond
    
    quartile_r_1 <- quantile(cv_r, 0.25)
	quartile_r_3 <- quantile(cv_r, 0.75)
    quartile_t_1 <- quantile(cv_t, 0.25)
	quartile_t_3 <- quantile(cv_t, 0.75)
	low_r_IQR <- quartile_r_1 - (1.5 * (quartile_r_3 - quartile_r_1))
    high_r_IQR <- quartile_r_3 + (1.5 * (quartile_r_3 - quartile_r_1))
    low_t_IQR <- quartile_t_1 - (1.5 * (quartile_t_3 - quartile_t_1))
    high_t_IQR <- quartile_t_3 + (1.5 * (quartile_t_3 - quartile_t_1))
    
    df <- data.frame(Median = c(median(cv_r), median(cv_t)),
                      Min = c(low_r_IQR, low_t_IQR),
                      Lower = c(quartile_r_1, quartile_t_1),
                      Upper = c(quartile_r_3, quartile_t_3),
                      Max = c(high_r_IQR, high_t_IQR),
                      Group = c(r_name, t_name))

    
    if (output$args$mean_correct) {
      
      cv_r_bg <- unlist(output$bg_metrics$means_r)
      cv_t_bg <- unlist(output$bg_metrics$means_t)
      
      r_bg_name <- paste(r_name, "(BG)")
      t_bg_name <- paste(t_name, "(BG)")
      
     
      quartile_rbg_1 <- quantile(cv_r_bg, 0.25)
	  quartile_rbg_3 <- quantile(cv_r_bg, 0.75)
      quartile_tbg_1 <- quantile(cv_t_bg, 0.25)
	  quartile_tbg_3 <- quantile(cv_t_bg, 0.75)

	  low_rbg_IQR <- quartile_rbg_1 - (1.5 * (quartile_rbg_3 - quartile_rbg_1))
      high_rbg_IQR <- quartile_rbg_3 + (1.5 * (quartile_rbg_3 - quartile_rbg_1))
      low_tbg_IQR <- quartile_tbg_1 - (1.5 * (quartile_tbg_3 - quartile_tbg_1))
      high_tbg_IQR <- quartile_tbg_3 + (1.5 * (quartile_tbg_3 - quartile_tbg_1))
    
      df <- rbind(df, data.frame(Median = c(median(cv_r_bg), median(cv_t_bg)),
                        Min = c(low_rbg_IQR, low_tbg_IQR),
                        Lower = c(quartile_rbg_1, quartile_tbg_1),
                        Upper = c(quartile_rbg_3, quartile_tbg_3),
                        Max = c(high_rbg_IQR, high_tbg_IQR),
                        Group = c(r_bg_name, t_bg_name)))
    }
    
    # Wes Anderson: Darjeeling
    # Blue Red Teal Orange
    palette <- c("#5BBCD6", "#FF0000", "#00A08A", "#F2AD00")
    
    df %>%
      ggplot(aes(x = Group, fill = Group)) +
      geom_boxplot(aes(ymin=Min, lower=Lower, middle=Median, upper=Upper, ymax=Max, group = Group), stat = "identity") +
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
