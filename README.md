# ConAn
Differential Network Connectivity Analysis

<img src="media/logo.png" height="100px" align="right"/>

[![](https://lifecycle.r-lib.org/articles/figures/lifecycle-maturing.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![](https://img.shields.io/github/last-commit/montilab/ConAn.svg)](https://github.com/montilab/ConAn/commits/master)

## Documentation

Please visit <https://montilab.github.io/ConAn>

## Requirements

We recommend the latest version of R (\>= 4.0.0) but **ConAn** currently
requires R (\>= 3.6.0) to be installed directly from Github.

### Install with devtools
```R
library(devtools)
devtools::install_github("montilab/ConAn")
```

```{r}
library(ConAn)
```

### Example
```R
simulations <- readRDS(file.path(system.file("extdata", package="ConAn"), "simulations.rds"))

eset <- simulations$eset
pathways <- simulations$pathways
names(pathways) <- paste("P", seq(length(pathways)), sep="")

output <- conan(eset=eset,
                mod_list=pathways,
                covariate="group",
                ctrl="G1",
                cond="G2",
                sim_type=c("bootstrap", "permutation")[[1]],
                iter=20,
                mean_correct=TRUE,
                cores=1,
                mdc_type=c("fraction", "difference")[[2]],
                plotting=TRUE,
                reporting=TRUE,
                report_path="report.Rmd")

```
