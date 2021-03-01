# ConAn
Differential Network Connectivity Analysis

### Install with devtools
```R
library(devtools)
devtools::install_github("anfederico/ConAn")
```

```{r}
library(ConAn)
```

### Run ConAn
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
                plotting=T,
                reporting=T,
                report_path="report.Rmd")

```

### Note: ConAn is under active development