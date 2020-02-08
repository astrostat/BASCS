mcmc.extended.full <- function(fix_runs,online_ordering,rjmcmc_run,w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num){
  
  # Standard MCMC updates
  for (t2 in 1:fix_runs){
    
    # Update photon allocations 
    probs <- matrix(NA,mix_num,obs_num)
    probs[1:k_curr,] <- t(matrix(unlist(lapply(1:k_curr,function(i) w[i]*psf_cpp(off.angle,ellip,slope,psf.norm,r0,spatial,mu_curr[,i])*(ewt[i]*dgamma(energy,eparas[i,3],eparas[i,3]/eparas[i,1])+(1-ewt[i])*dgamma(energy,eparas[i,4],eparas[i,4]/eparas[i,2])))),ncol=k_curr))
    probs[mix_num,energy <= max_back_energy] <- w[mix_num]*(1/img_area)*(1/max_back_energy)
    probs[mix_num,energy > max_back_energy] <- 0
    allocate_curr <- t(matrix(unlist(lapply(1:obs_num,function(i) rmultinom(1, 1, probs[,i]))),ncol=obs_num))  # Don't need to normalize as rmultinom function does it automatically

    # Counts vector
    mix_num <- k_curr+1
    count_vector <- matrix(NA,mix_num,1)
    count_vector[1:mix_num] <- apply(allocate_curr[,1:mix_num],2,sum)
    
    # Update positions
    mu_prop <- mu_curr 
    for (i in 1:k_curr){      
      index <- allocate_curr[,i]==1
      if (count_vector[i]>0){
        # Adaptive version (eventually ended to ensure convegence)
        if (rjmcmc_run < adapt_end){
          mu_prop[,i] <- rnorm(2,mu_curr[,i],mu_adapt_prop_sd/sqrt(count_vector[i]))
        }
        # Non-adaptive version
        if (rjmcmc_run >= adapt_end){
          mu_prop[,i] <- rnorm(2,mu_curr[,i],mu_fixed_prop_sd)
        }
        logr <- sum(log(psf_cpp(off.angle,ellip,slope,psf.norm,r0,spatial[index,],mu_prop[,i])))-sum(log(psf_cpp(off.angle,ellip,slope,psf.norm,r0,spatial[index,],mu_curr[,i])))
        u <- runif(1,0,1)
        if(is.na(logr)==0){
          if (log(u) < logr){
            mu_curr[,i] <- mu_prop[,i]
          }
        }
      }
      # Have to make sure that sources without phoons assigned move (doesn't effect likelihood) 
      if (count_vector[i]==0){
        if (rjmcmc_run < adapt_end){
          mu_curr[,i] <- c(runif(1,xlow,xup),runif(1,ylow,yup))
        } else {
          mu_curr[,i] <- rnorm(2,mu_curr[,i],mu_fixed_prop_sd)
        }
      }
    }
    
    # Order parameters by suspected sources intensities (associated with particular positions)
    if (k_curr > 1 & online_ordering =="reference"){
      to_order <- min(no_guess,k_curr)
      next_index <- which.min(apply((mu_curr-mu_guess[,1])^2,2,sum))
      next_index_store <- next_index
      if (to_order > 1){
        for (i in 2:to_order){
          next_order <- order(apply((mu_curr-mu_guess[,i])^2,2,sum))
          next_index_all <- setdiff(next_order,next_index_store)
          next_index_store <- c(next_index_store,next_index_all[1])
        }
      }
      indexmu <- c(next_index_store,setdiff(1:k_curr,next_index_store))
      mu_curr <- mu_curr[,indexmu]
      count_vector[1:k_curr] <- count_vector[indexmu]
      allocate_curr[,1:k_curr] <- allocate_curr[,indexmu]
      eparas <- eparas[indexmu,] 
      ewt <- ewt[indexmu]
    }
    
    # Update weights
    alpha <- rep(wprior,mix_num)
    w <- rdirichlet(1,alpha+count_vector)
    
    # Update spectral parameters (extended full model) 
    for (i in 1:k_curr){
      index <- which(allocate_curr[,i] == 1)
      cspatial <- energy[index]
      gm1curr <- eparas[i,1]
      gm2curr <- eparas[i,2]
      ga1curr <- eparas[i,3]
      ga2curr <- eparas[i,4]
      ewtcurr <- ewt[i]
      ewtprop <- inv.logit(rnorm(1,logit(ewtcurr),specwt_sd))
      gm1prop <- rnorm(1,gm1curr,i*specm_sd)
      ga1prop <- rnorm(1,ga1curr,speca_sd)
      if ((gm1prop > emean.min) & (gm1prop < emean.max) & (ga1prop > 0)){
        logr <- spectral_post(gm1prop,ga1prop,gm2curr,ga2curr,ewtprop,cspatial) - spectral_post(gm1curr,ga1curr,gm2curr,ga2curr,ewtcurr,cspatial)
        u <- runif(1,0,1)
        if (log(u) < logr){
          eparas[i,c(1,3)] <- c(gm1prop,ga1prop)
          ewt[i] <- ewtprop
          gm1curr <- eparas[i,1]
          ga1curr <- eparas[i,3]
          ewtcurr <- ewt[i]
        }
      }
      gm2prop <- rnorm(1,gm2curr,i*specm_sd)
      ga2prop <- rnorm(1,ga2curr,speca_sd)
      if ((gm2prop > emean.min) & (gm2prop < emean.max) & (ga2prop > 0)){
        logr <- spectral_post(gm1curr,ga1curr,gm2prop,ga2prop,ewtcurr,cspatial) - spectral_post(gm1curr,ga1curr,gm2curr,ga2curr,ewtcurr,cspatial)
        u <- runif(1,0,1)
        if (log(u) < logr){
          eparas[i,c(2,4)] <- c(gm2prop,ga2prop)
        }
      }
    }
    # Order Gammas for identifiability and combine conditions
    if (k_curr > 1){
      epara_order <- apply(eparas[,c(1,2)],1,order)
      epara_order_all <- rbind(epara_order,epara_order+2)
      for (i in 1:k_curr){
        eparas[i,] <- eparas[i,epara_order_all[,i]] 
      }
      for (i in 1:k_curr){
        if (sum(epara_order[,i] != c(1,2))>0){
          ewt[i] <- 1- ewt[i]
        }
      }
    }

    
  }
  
  # Output parameters and log-posterior
  value <- list(c(k_curr,c(mu_curr),c(w),c(eparas),c(ewt)),allocate_curr,extended_full_log_posterior(mu_curr,allocate_curr,w,k_curr,eparas,ewt,alpha))
  return(value)
}
  