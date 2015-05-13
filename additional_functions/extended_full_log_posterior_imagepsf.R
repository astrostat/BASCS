extended_full_log_posterior <- function(mu_curr,allocate_curr,w,k_curr,eparas,ewt,alpha){
  mix_num <- k_curr+1
  loglike <- 0
  like_obs <- matrix(NA,obs_num,mix_num)
  for (j in 1:k_curr){
    like_obs[,j] <- w[j]*psf(spatial,mu_curr[,j])*(ewt[j]*dgamma(energy,eparas[j,3],eparas[j,3]/eparas[j,1])+(1-ewt[j])*dgamma(energy,eparas[j,4],eparas[j,4]/eparas[j,2]))
  }
  like_obs[,mix_num] <- (w[mix_num]/img_area)*(1/max_back_energy)
  loglike <- sum(log(apply(like_obs,1,sum)))
  components <- matrix(NA,5,1)
  components[1] <- loglike # log likelihood
  components[2] <- -k_curr*log(img_area)   # mu_prior
  components[3] <- log(ddirichlet(w,alpha))   # w_prior
  components[4] <- log(dpois(k_curr,theta))   # k_prior 
  components[5] <- sum(log(dgamma(eparas[,c(3,4)],ashape,arate)))-2*k_curr*log(emean.range) + sum(log(dbeta(ewt,ewt_ab,ewt_ab))) # eparas_prior
  value <- c(sum(components),loglike)
  return(value)
}