dis_lik <- function(abund = NULL,
                    bands = NULL,
                    max_dist = NULL,
                    max_intervals = NULL,
                    tau = NULL)
{
  loglik <- vector(mode = "double", length = length(bands))
  
  for (i in 1:length(bands))
  {
    Pi <- vector(mode = "double", length = bands[i])
    for (k in 1:(bands[i] - 1))
    {
      if (k > 1)
      {
        Pi[k] <- ((1 - exp(-(max_dist[i,k]^2 / (tau[i])^2))) - 
                    (1 - exp(-(max_dist[i,k - 1]^2 / (tau[i])^2)))) / 
          (1 - exp(-(max_dist[i,bands[i]]^2 / (tau[i])^2)))
      }else{
        Pi[k] = (1 - exp(-(max_dist[i,k]^2 / (tau[i])^2))) /
          (1 - exp(-(max_dist[i,bands[i]]^2 / (tau[i])^2)))
      }
    }
    
    Pi[bands[i]] = 1 - sum(Pi)
    
    loglik[i] <- dmultinom(x = abund[i,1:bands[i]], prob = Pi)
  }
  
  return(loglik)
}