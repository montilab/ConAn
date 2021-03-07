#' @import magrittr
#' @import hypeR
#' @import ggpubr
#' @import kableExtra
#' @keywords internal
report <- function(output) {
  
  genesets <- msigdb_gsets("Homo sapiens", "C2", "CP:KEGG", clean=TRUE)
  mhyp <- hypeR(mod_list, genesets, test="hypergeometric", background=30000)
  REACTOME <- msigdb_gsets(species="Homo sapiens", category="C2", subcategory="CP:REACTOME")
  mhyp_R <- hypeR(mod_list, REACTOME, test="hypergeometric", background=30000)
  
  r_name <- output$args$ctrl
  t_name <- output$args$cond

  mdc <- as.data.frame(output$data$mod_n)

  mdc <- cbind(mdc, as.data.frame(output$stat[c("mods_mc_r",
                                                "mods_mc_t",
                                                "mods_mdc_adj")]))
  mdc <- cbind(mdc, output$significance)

  if (output$args$sim_type == "bootstrap") {
      colnames(mdc) <- c("Gene Size",
                         paste(r_name, "Connectivity"),
                         paste(t_name, "Connectivity"),
                         "MDC",
                         "Stat",
                         "Stdev",
                         "P-Value",
                         "FDR")
  } else {
      colnames(mdc) <- c("Gene Size",
                         paste(r_name, "Connectivity"),
                         paste(t_name, "Connectivity"),
                         "MDC",
                         "P-Value",
                         "FDR")
  }
  
  sigmdc <- mdc %>% 
    filter(FDR <= 0.05)

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
```{r setup, echo=FALSE, message=FALSE}
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
ggarrange(p1, p2, ncol=2, widths=c(0.4, 0.6))
```
"

  rmd_results <- "
# Results
```{r}
mdc
```
"

  rmd_hyp <- "
### hypeR Enrichment Analysis and Significantly Differential Modules
```{r {2}, fig.width=15, fig.align='center'}
p2 <- hyp_dots(mhyp, merge=TRUE, fdr=0.05, title='KEGG')
p3 <- hyp_dots(mhyp_R, merge=TRUE, fdr=0.05, title='Reactome')
ggarrange(p2, p3, ncol=2, widths=c(0.5, 0.5))
```
"

  rmd_sigresults <- "
# Results
```{r}
sigmdc
```
"


  file_path <- output$args$report_path

  write(rmd_config, file=file_path, append=FALSE)
  write(rmd_arguments, file=file_path, append=TRUE)
  write(rmd_knitr, file=file_path, append=TRUE)
  if (output$args$plotting) {
      write(rmd_tabset, file=file_path, append=TRUE)

      for (tab in output$data$mod_names) {
          rmd_tab %>%
          format_str(tab) %>%
          write(file=file_path, append=TRUE)
      }
  }
  write(rmd_results, file=file_path, append=TRUE)
  write(rmd_hyp, file=file_path, append=TRUE)
  write(rmd_sigresults, file=file_path, append=TRUE)
  
  rmarkdown::render(input=file_path,
                    output_format="html_document",
                    output_file=paste(file_path, "html", sep="."))
}
  

