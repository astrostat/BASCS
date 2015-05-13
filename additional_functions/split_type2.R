# Extended full model split move type 2 
split_type2 <- function(w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num,count_zero){
  
  # Initializae acceptance probability
  A <- 0
  
  # Intialize failure diagnostic
  fail_reason <- 1
  
  # Draw variables for dimension matching
  s1p <- sample(c(-1,1),1)
  s2p <- sample(c(-1,1),1)
  u2 <- s1p*rbeta(1,2,2)
  u3 <- s2p*rbeta(1,2,2)
  # Choose component to split
  pick <- sample(1:k_curr,1)
  # Performance check
  if (sum(allocate_curr[,pick]==1)==0){
    count_zero <- count_zero+1
  }
  # Calculate proposal parameters
  w.com <- w[pick]
  u1 <- ewt[pick]
  w1 <- w.com*u1
  w2 <- w.com*(1-u1)
  mu.com <- mu_curr[,pick]
  mu1 <- c(mu.com[1] - u2*sqrt(sigma[1,1]*(1-u1)/u1),mu.com[2] - u3*sqrt(sigma[1,1]*(1-u1)/u1))
  mu2 <- c(mu.com[1] + u2*sqrt(sigma[1,1]*u1/(1-u1)),mu.com[2] + u3*sqrt(sigma[1,1]*u1/(1-u1))) 
  
  # Check if mu2 is closest source to mu1, and if the positions are within the range of the data
  if ((sum((mu1-mu2)^2) < min(apply((mu1-mu_curr[,-pick])^2,2,sum))) & (max(mu1[1],mu2[1]) < xup) & (min(mu1[1],mu2[1]) > xlow) & (max(mu1[2],mu2[2]) < yup) & (min(mu1[2],mu2[2]) > ylow)){
    fail_reason <- NA
    
    ##############################
    # Energy parameters
    ewt.com <- ewt[pick]
    # New spectral weights (heavily favor first component)
    ewt1 <- rbeta(1,prop_ewta,prop_ewtb)
    ewt2 <- rbeta(1,prop_ewta,prop_ewtb)
    
    # Origianl parameters
    gam.com.mean1 <- eparas[pick,1]
    gam.com.shape1 <- eparas[pick,3]
    gam.com.mean2 <- eparas[pick,2]
    gam.com.shape2 <- eparas[pick,4]
    
    # First source parameters
    gam.mean11 <- gam.com.mean1
    gam.shape11 <- gam.com.shape1
    v2 <- runif(1,0,1)
    gam.mean12 <- gam.com.mean1 + v2*(emean.max-gam.com.mean1)
    v4 <- rbeta(1,1,s)
    gam.shape12 <- a_scale*v4
    
    # Second source parameters
    gam.mean21 <- gam.com.mean2
    gam.shape21 <- gam.com.shape2
    v3 <- runif(1,0,1)
    gam.mean22 <- gam.com.mean2 + v3*(emean.max-gam.com.mean2)
    v5 <- rbeta(1,1,s)
    gam.shape22 <- a_scale*v5
    ##############################
    
    
    # Define proposal parameter vectors
    # Weights
    w.new <- matrix(NA,k_curr+2,1)
    w.new[setdiff(1:(k_curr+2),pick:(pick+1))] <- w[-pick]
    w.new[pick:(pick+1)] <- c(w1,w2)
    # Positions
    mu.new <- matrix(NA,2,k_curr+1)
    mu.new[,setdiff(1:(k_curr+1),pick:(pick+1))] <- mu_curr[,-pick]
    mu.new[,pick] <- mu1
    mu.new[,(pick+1)] <- mu2
    # Eparas
    ewt.new <-  matrix(NA,k_curr+1,1)
    ewt.new[setdiff(1:(k_curr+1),pick:(pick+1))] <- ewt[-pick]
    ewt.new[pick] <- ewt1
    ewt.new[pick+1] <- ewt2
    eps.new <- matrix(NA,k_curr+1,4)
    eps.new[setdiff(1:(k_curr+1),pick:(pick+1)),] <- eparas[-pick,]
    eps.new[pick,] <- c(gam.mean11,gam.mean12,gam.shape11,gam.shape12)
    eps.new[(pick+1),] <- c(gam.mean21,gam.mean22,gam.shape21,gam.shape22)
    # Specify proposed allcoation vector (randomly allocate points from split component to two new 
    # components based on parameters)
    new.allocate <- matrix(0,obs_num,k_curr+2)
    new.allocate[,-pick] <- allocate_curr
    # Allocation and allocation probabilities
    if (sum(allocate_curr[,pick]==1)>0){
      indcom <- which(allocate_curr[,pick]==1)
      probs.new <- matrix(NA,1,length(indcom))
      p1 <- w1*psf_cpp(off.angle,ellip,slope,psf.norm,r0,spatial[indcom,],mu1)*(ewt1*dgamma(energy[indcom],gam.shape11,gam.shape11/gam.mean11)+(1-ewt1)*dgamma(energy[indcom],gam.shape12,gam.shape12/gam.mean12))
      p2 <- w2*psf_cpp(off.angle,ellip,slope,psf.norm,r0,spatial[indcom,],mu2)*(ewt2*dgamma(energy[indcom],gam.shape21,gam.shape21/gam.mean21)+(1-ewt2)*dgamma(energy[indcom],gam.shape22,gam.shape22/gam.mean22))
      probs.new <- p1/(p1+p2)
      new.allocate[indcom,pick] <- rbinom(length(indcom),1,probs.new)
      new.allocate[,pick+1] <- allocate_curr[,pick]-new.allocate[,pick]
    }
    # Acceptance probability
    A <- exp(log_accept_split_extended_full_type2(k_curr,mu_curr,mu.new,allocate_curr,new.allocate,pick,pick+1,eparas,eps.new,ewt,ewt.new,pick,w,w.new,probs.new,c(ewt.com,u2,u3,ewt1,ewt2,v4,v5)))
    
    utest <- runif(1,0,1)
    if (is.na(A)==0){
      if (utest < A){
        w <- w.new
        allocate_curr <- new.allocate
        mu_curr <- mu.new
        eparas <- eps.new
        ewt <- ewt.new
        k_curr <- k_curr+1
        mix_num <- k_curr + 1
      }
    }
  }
  
  value <- list(w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num,fail_reason,mu1,mu2,A)
  return(value)
}