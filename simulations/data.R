genes <- paste("G", seq(300), sep="")

# Shuffle pathways
i <- 1
pathways <- list()
set.seed(1)
while (length(genes) > 45) {
    pathways[[i]] <- sample(genes, sample(seq(15,45), 1), replace=F)
    genes <- setdiff(genes, unlist(pathways))
    i <- i+1
}
pathways[[i]] <- genes
table(unlist(lapply(pathways, length)))

# Reset
genes <- paste("G", seq(300), sep="")
table(sort(unlist(pathways)) == sort(genes))

# True Networks
adj.1 <- matrix(0, nrow=length(genes), ncol=length(genes), dimnames=list(genes, genes))
adj.2 <- matrix(0, nrow=length(genes), ncol=length(genes), dimnames=list(genes, genes))

true.connectivity <- list(c(0.1, 0.95),
                          c(0.99, 0.01),
                          c(0.4, 0.75),
                          c(0.5, 0.5),
                          c(0.5, 0.5),
                          c(0.5, 0.5),
                          c(0.5, 0.5),
                          c(0.5, 0.5),
                          c(0.5, 0.5),
                          c(0.5, 0.5))

for (i in seq(length(pathways))) {
    p <- pathways[[i]]
    adj.1[p, p] <- rbinom(length(p)^2, 9, true.connectivity[[i]][[1]])/10
    adj.2[p, p] <- rbinom(length(p)^2, 9, true.connectivity[[i]][[2]])/10
}

# Make symmetric
diag(adj.1) <- 1
diag(adj.2) <- 1
adj.1[lower.tri(adj.1)] = t(adj.1)[lower.tri(adj.1)]
adj.2[lower.tri(adj.2)] = t(adj.2)[lower.tri(adj.2)]

sh <- lapply(pathways, function(p) {
    print(mean(adj.1[p, p]) / mean(adj.2[p, p])) 
})

# Make positive definite
library(lqmm)
adj.1.pd <- make.positive.definite(adj.1)
adj.2.pd <- make.positive.definite(adj.2)

sh <- lapply(pathways, function(p) {
    print(mean(adj.1.pd[p, p]) / mean(adj.2.pd[p, p])) 
})

# Generate multivariate guassian networks
library(MASS)
n <- 20
mvg.1 <- mvrnorm(n, mu=rep(0, length(genes)), Sigma=adj.1.pd)
mvg.2 <- mvrnorm(n, mu=rep(0, length(genes)), Sigma=adj.2.pd)

# Expression Set
library(Biobase)
library(magrittr)
fdat <- data.frame(symbol=genes) %>% set_rownames(genes)
edat <- cbind(t(mvg.1), t(mvg.2)) %>% set_rownames(genes) %>% set_colnames(seq(ncol(.)))
pdat <- data.frame(group=c(rep("G1", n), rep("G2", n)))
eset <- ExpressionSet(edat,
                      AnnotatedDataFrame(pdat),
                      AnnotatedDataFrame(fdat))

sh <- lapply(pathways, function(p) {
    edat.1 <- exprs(eset[p, eset$group == "G1"])
    edat.2 <- exprs(eset[p, eset$group == "G2"])
    print(mean(cor(t(edat.1))) / mean(cor(t(edat.2)))) 
})

simulations <- list(eset=eset,
                    pathways=pathways)

saveRDS(simulations, "simulations.rds")
