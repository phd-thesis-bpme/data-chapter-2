rem_lik <- function(abund = NULL,
                      bands = NULL,
                      max_time = NULL,
                      max_intervals = NULL,
                      phi = NULL)
{
  lp <- vector(mode = "double", length = length(bands))
  
  for (i in 1:length(bands))
  {
    Pi <- vector(mode = "double", length = bands[i])
    for (j in 2:bands[i])
    {
      Pi[j] <- (exp(-max_time[i,j-1] * phi[i]) -
                    exp(-max_time[i,j] * phi[i])) /
        (1 - exp(-max_time[i,bands[i]] * phi[i]))
    }
    Pi[1] = 1 - sum(Pi)
    
    lp[i] <- dmultinom(x = abund[i,1:bands[i]], prob = Pi)
  }
  
  return(lp)
}