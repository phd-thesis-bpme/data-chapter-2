generate_phylo_corr <- function(phylo_tree = NULL, drops = NULL)
{
  # This code adopted from Solymos et al. 2018
  CORR <- TRUE
  vv <- list()
  vv[[1]] <- vcv(drop.tip(phylo_tree[[1]], tip = drops), corr = CORR)
  for (i in 2:length(phylo_tree))
  {
    v <- ape::vcv(drop.tip(phylo_tree[[i]], tip = drops), corr = CORR)
    v <- v[rownames(vv[[1]]), colnames(vv[[1]])]
    vv[[i]] <- v
  }
  vvv <- v
  for (i in 1:length(v))
  {
    vvv[i] <- mean(sapply(vv, function(z) z[i]))
  }
  
  return(as.matrix(Matrix::nearPD(vvv, corr = CORR)$mat))
}
