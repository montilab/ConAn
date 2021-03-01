simulations <- readRDS("simulations/simulations.rds")

eset <- simulations$eset
pathways <- simulations$pathways
names(pathways) <- paste("P", seq(length(pathways)), sep="")

# Run Differential Connectivity Analysis
r1 <- conan(eset=eset,
                mod_list=pathways,
                covariate="group",
                ctrl="G1",
                cond="G2",
                sim_type=c("bootstrap", "permutation")[[1]],
                iter=100,
                mean_correct=F,
                cores=3,
                mdc_type=c("frac", "diff")[[1]])

r2 <- conan(eset=eset,
                mod_list=pathways,
                covariate="group",
                ctrl="G1",
                cond="G2",
                sim_type=c("bootstrap", "permutation")[[2]],
                iter=100,
                mean_correct=F,
                cores=3,
                mdc_type=c("frac", "diff")[[1]])

r3 <- conan(eset=eset,
                mod_list=pathways,
                covariate="group",
                ctrl="G1",
                cond="G2",
                sim_type=c("bootstrap", "permutation")[[1]],
                iter=100,
                mean_correct=F,
                cores=3,
                mdc_type=c("frac", "diff")[[2]])

r4 <- conan(eset=eset,
                mod_list=pathways,
                covariate="group",
                ctrl="G1",
                cond="G2",
                sim_type=c("bootstrap", "permutation")[[2]],
                iter=100,
                mean_correct=F,
                cores=3,
                mdc_type=c("frac", "diff")[[2]])

r5 <- conan(eset=eset,
                mod_list=pathways,
                covariate="group",
                ctrl="G1",
                cond="G2",
                sim_type=c("bootstrap", "permutation")[[1]],
                iter=100,
                mean_correct=T,
                cores=3,
                mdc_type=c("frac", "diff")[[1]])

r6 <- conan(eset=eset,
                mod_list=pathways,
                covariate="group",
                ctrl="G1",
                cond="G2",
                sim_type=c("bootstrap", "permutation")[[2]],
                iter=100,
                mean_correct=T,
                cores=3,
                mdc_type=c("frac", "diff")[[1]])

r7 <- conan(eset=eset,
                mod_list=pathways,
                covariate="group",
                ctrl="G1",
                cond="G2",
                sim_type=c("bootstrap", "permutation")[[1]],
                iter=100,
                mean_correct=T,
                cores=3,
                mdc_type=c("frac", "diff")[[2]])

r8 <- conan(eset=eset,
                mod_list=pathways,
                covariate="group",
                ctrl="G1",
                cond="G2",
                sim_type=c("bootstrap", "permutation")[[2]],
                iter=100,
                mean_correct=T,
                cores=3,
                mdc_type=c("frac", "diff")[[2]])

data <- list(r1,r2,r3,r4,r5,r6,r7,r8)

ground.truth <- readRDS("dev/data/Ground-Truth.rds")

for (i in seq(length(data))) {
    print(table(data[[i]]$stat$mods_mdc_adj == ground.truth[[i]]$stat$mods_mdc_adj))
}
