simulations <- readRDS("simulations/simulations.rds")

eset <- simulations$eset
pathways <- simulations$pathways
names(pathways) <- paste("P", seq(length(pathways)), sep="")

#simulations <- readRDS(file.path(system.file("extdata", package="ConAn"), "simulations.rds"))

#eset <- simulations$eset
#pathways <- simulations$pathways
#names(pathways) <- paste("P", seq(length(pathways)), sep="")

exprs(eset)[1,] <- rep(1, 40)
exprs(eset)[5,] <- rep(0, 40)
apply(exprs(eset), 1, var)

pathways$P10 <- c(pathways$P10, c("Z1", "Z2", "Z3", "Z4"))

# Run Differential Connectivity Analysis
output <- conan(eset=eset,
                mod_list=pathways,
                covariate="group",
                ctrl="G1",
                cond="G2",
                sim_type=c("bootstrap", "permutation")[[2]],
                iter=20,
                mean_correct=TRUE,
                cores=1,
                mdc_type=c("fraction", "difference")[[2]],
                plotting=T,
                reporting=T,
                report_path="report.Rmd")


#P1       P2       P3       P4       P5       P6       P7       P8       P9      P10 
#0.64963 -0.66625  0.20774 -0.03062  0.09110 -0.13047  0.11436  0.11167  0.00582 -0.07804 

output$significance


library(wesanderson)

wes_palette("Darjeeling1")


# Wes Anderson: Darjeeling
# Red Blue Orange Teal
palette <- c("#FF0000", "#5BBCD6", "#F2AD00", "#00A08A")
