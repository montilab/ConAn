# ConAn
Differential Network Connectivity Analysis

### Install with devtools
```R
library(devtools)
devtools::install_github("montilab/ConAn", auth=<Personal Access Token>)

# For help with generating Personal Access Tokens
# https://help.github.com/articles/creating-a-personal-access-token-for-the-command-line/
```

### Run ConAn
```R
library(ConAn)

# Load test data
eset <- ConAn::eset
mod_list <- ConAn::mod_list

# Run ConAn
output <- ConAn::diff_conn(eset,
                           covariate = "subtype",
                           ctrl = "LumA",
                           cond = "Basal",
                           sim_type = "bootstrap",
                           iter = 5,
                           mean_correct = T,
                           cores = 1,
                           use_gpu = F,
                           mdc_type = "frac",
                           plotting = T)

# Plot data
output$plots
```

### Note: ConAn is under active development