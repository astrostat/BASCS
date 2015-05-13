# Extended full model birth move
birth_move_extended_full <- function(w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num){
  
  # Count the number of empty sources
  count_short_new <- ifelse(k_curr > 1,apply(allocate_curr[,1:k_curr],2,sum),sum(allocate_curr[,1]))
  empty <- sum(count_short_new == 0)
  
  # Proposal - except for w, draw all parameters from their priors
  wbirth <- rbeta(1,abirth,bbirth)
  munew <- matrix(NA,2,k_curr+1)
  munew[,1:k_curr] <- mu_curr
  additional_source <- c(runif(1,xlow,xup),runif(1,ylow,yup))
  munew[,k_curr+1] <- additional_source
  wnew <- matrix(NA,1,k_curr+2)
  wnew[-(k_curr+1)] <- w*(1-wbirth)
  wnew[k_curr+1] <- wbirth
  eps_new <- matrix(NA,k_curr+1,4)
  eps_new[1:k_curr,] <- eparas
  enew <- c(runif(2,emean.min,emean.max),rgamma(2,ashape,arate))
  # Need to order gamma means 
  # Note ewt_new[k_curr+1] must be with smaller component as that is how it is in the prior 
  # and we are assuming the prior cancels with proposal
  if (enew[1] < enew[2]){
    eps_new[k_curr+1,] <- enew
  } else {
    eps_new[k_curr+1,] <- enew[c(2,1,4,3)]
  }
  ewt_new <- matrix(NA,k_curr+1,1)
  ewt_new[1:k_curr] <- ewt
  ewt_new[k_curr+1] <- rbeta(1,ewt_ab,ewt_ab)
  
  # Acceptance probability
  A <- exp(log_accept_birth(k_curr,empty,wbirth))
  
  utest <- runif(1,0,1)
  if (is.na(A)==0){
    if (utest < A){
      mu_curr <- munew
      eparas <- eps_new
      ewt <- ewt_new
      k_curr <- k_curr+1
      mix_num <- k_curr+1
      w <- wnew
      new_allocate <- matrix(0,obs_num,mix_num)
      new_allocate[,-k_curr] <- allocate_curr
      allocate_curr <- new_allocate
    }
  }
  
  value <- list(w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num,additional_source,A)
  return(value)
}