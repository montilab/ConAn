#' @title Differential Connectivity Report
#' @description Generate an rmarkdown report of the results returned by conan
#' @param output The object returned by a call to conan
#' @param FDR_thresh Max fdr to report as significant
#' @param mods The modules to report on (default is all)
#' @param file_path Where to save the rmarkdown (an '.html' will be generated as well)
#' @param digits Number of significant digits to display in the tables
#'
#' @import magrittr
#' @import ggpubr
#' @import kableExtra
#' @import DT
#'
#' @export
report <- function(
  output, 
  FDR_thresh, 
  mods=output$data$mod_list, 
  file_path,
  digits=2,
  annotate_mods=FALSE
) 
{
  ## check inputs
  if (FDR_thresh<=0 || FDR_thresh>=1.0) stop( "FDR_thresh must be in (0,1)")
  
  r_name <- output$args$ctrl
  t_name <- output$args$cond

  ## Module Differential Connectivity table
  mdc <- as.data.frame(output$data$mod_n)
  mdc <- cbind(mdc, as.data.frame(output$stat[c("mods_mc_r",
                                                "mods_mc_t",
                                                "mods_mdc_adj")]))
  mdc <- cbind(mdc, output$significance) %>%
    dplyr::mutate(across(where(is.numeric),signif,digits=digits))

  if (output$args$sim_type == "bootstrap") {
      colnames(mdc) <- c("Module Size",
                         paste(r_name, "\nConnectivity"),
                         paste(t_name, "\nConnectivity"),
                         "MDC",
                         "Stat",
                         "Stdev",
                         "P-Value",
                         "FDR")
  } else {
      colnames(mdc) <- c("Module Size",
                         paste(r_name, "\nConnectivity"),
                         paste(t_name, "\nConnectivity"),
                         "MDC",
                         "P-Value",
                         "FDR")
  }
  ## MDC table restricted to significant modules
  sigmdc <- mdc %>% 
    dplyr::filter(FDR <= FDR_thresh) %>%
    dplyr::mutate(across(where(is.numeric),signif,digits=digits))

  ## Annotate modules with fData
  if ( annotate_mods && ncol(Biobase::fData(output$data$eset))>0 ) {
    mods <- lapply(mods, function(M) Biobase::fData(output$data$eset)[M,])
  }
  ## Start markdown generation
  rmd_config <- "---
title: 'Differential Connectivity Analysis Report'
date: 'Report Created: `r Sys.Date()`'
output:
  html_document:
    theme: simplex
    toc: false
    df_print: paged
---
"

  rmd_knitr <- "
```{r conan.setup, echo=FALSE, message=FALSE}
knitr::opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE)
library(kableExtra)
options(scipen=1, digits=3)
```
"

  rmd_arguments <- "
# Parameters
**Covariate**: `r output$args$covariate`  
**Reference Group**: `r output$args$ctrl`  
**Test Group**: `r output$args$cond`  
**Simulation Type**: `r output$args$sim_type`  
**Iterations**: `r output$args$iter`  
**Cores**: `r output$args$cores`  
**Mean Background Correction**: `r output$args$mean_correct`  
**Background Sample Size**: `r output$args$N_genes`  
**Background Iterations**: `r output$args$iter_bg`  
**Differential Connectivity Type**: `r output$args$mdc_type`  
**Plotting**: `r output$args$plotting`  
**Report Directory**: `r output$args$report_path`  

***

# Background Statistics
**`r output$args$ctrl` Background**: `r output$bg$mc_r_bg`  
**`r output$args$cond` Background**: `r output$bg$mc_t_bg`  

***
"
  rmd_tabset <- "
# Visualizations
##  {.tabset .tabset-fade}
  "

  rmd_tab <- "
### {1} 
```{r {1}, fig.width=9, fig.align='center'}
p1 <- output$plots$connectivity[['{1}']]
p2 <- output$plots$permutations[['{1}']]
ggpubr::ggarrange(p1, p2, ncol=2, widths=c(0.4, 0.6))
```
"

  mod_tabset <- "
# Modules
##  {.tabset .tabset-fade}
  "
  
  mod_tab <- "
### {1}
```{r {1}_mod}
 DT::datatable(cbind(genes=mods[['{1}']]), options = list(scrollX = TRUE, paging=TRUE))
```
"

  rmd_results <- "
# Results
```{r}
 DT::datatable(mdc, options = list(scrollX = TRUE, paging=TRUE))
```
"

  rmd_sigresults <- "
# Significant Results
```{r}
 DT::datatable(sigmdc, options = list(scrollX = TRUE, paging=TRUE))
```
"

  write(rmd_config, file=file_path, append=FALSE)
  write(rmd_knitr, file=file_path, append=TRUE)
  write(rmd_arguments, file=file_path, append=TRUE)
  write(rmd_results, file=file_path, append=TRUE)
  write(rmd_sigresults, file=file_path, append=TRUE)
  
  if (output$args$plotting) {
      write(rmd_tabset, file=file_path, append=TRUE)

      for (tab in output$data$mod_names) {
          rmd_tab %>%
          format_str(tab) %>%
          write(file=file_path, append=TRUE)
      }
  }
  write(mod_tabset, file=file_path, append=TRUE)
  
  for (tab in output$data$mod_names) {
    mod_tab %>%
      format_str(tab) %>%
      write(file=file_path, append=TRUE)
  }
  html_file <- {
    if (stringr::str_ends(file_path,"Rmd"))
      stringr::str_replace(file_path,"Rmd","html")
    else
      paste(file_path, "html", sep=".")
  }
  rmarkdown::render(input=file_path,
                    output_format="html_document",
                    output_file=html_file)
}
