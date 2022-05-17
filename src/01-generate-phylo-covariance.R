####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 01-generate-phylo-covariance.R
# Created March 2022
# Last Updated April 2022

####### Import Libraries and External Files #######

library(ape)
library(Matrix)
library(ggplot2)

####### Read Data #################################

turdidae <- ape::read.nexus(file = "data/turdidae.nex")
binomial <- read.csv("output/binomial_names.csv")

####### Main Code #################################

# Consensus tree
tree_consensus <- ape::consensus(turdidae)

# This code adopted from Solymos et al. 2018
CORR <- TRUE
vv <- list()
vv[[1]] <- vcv(turdidae[[1]], corr = CORR)
for (i in 2:length(turdidae))
{
  v <- ape::vcv(turdidae[[i]], corr = CORR)
  v <- v[rownames(vv[[1]]), colnames(vv[[1]])]
  vv[[i]] <- v
}
vvv <- v
for (i in 1:length(v))
{
  vvv[i] <- mean(sapply(vv, function(z) z[i]))
}

corr_matrix <- as.matrix(Matrix::nearPD(vvv, corr = CORR)$mat)

####### Output ####################################

save(corr_matrix, file = "output/turdidae_corr_matrix.rda")

png(filename = "plots/turdidae_phylo_corr.png",
    width = 6, height = 6, units = "in", res = 300)
data <- cbind(expand.grid(dimnames(corr_matrix)), value = as.vector(corr_matrix))
print(ggplot(data = data, aes(x = Var1, y = Var2, fill = value)) +
        geom_tile())
dev.off()

png(filename = "plots/turdidae_phylo.png",
    width = 6, height = 6, units = "in", res = 300)
plot(tree_consensus)
dev.off()
