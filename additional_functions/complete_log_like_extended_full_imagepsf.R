# COMPLETE DATA spatial log-likelihood (conditional on missing spatial (and n))
complete_log_like_extended_full <- function(mu,allocate,eps,eweight,weight,k){
  loglike <- 0
  for (j in 1:k){
    index <- allocate[,j] == 1
    if (sum(index)>0){
      loglike <- loglike + sum(log(psf(spatial[index,],mu[,j]))) + sum(log(eweight[j]*dgamma(energy[index],eps[j,3],eps[j,3]/eps[j,1])+(1-eweight[j])*dgamma(energy[index],eps[j,4],eps[j,4]/eps[j,2])))
    }
  }
  loglike <- loglike + sum(allocate[,k+1])*(log(1/img_area) + log(1/max_back_energy))
  return(loglike)
}