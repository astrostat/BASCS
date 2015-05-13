# Extended full model split move type 1 
split_type1 <- function(w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num,count_zero){
  
  # Initializae acceptance probability
  A <- 0
  
  # Intialize failure diagnostic
  fail_reason <- 1
  
  # Draw variables for dimension matching
  s1p <- sample(c(-1,1),1)
  s2p <- sample(c(-1,1),1)
  u1 <- rbeta(1,u1a,u1b)
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
  w1 <- w.com*u1
  w2 <- w.com*(1-u1)
  mu.com <- mu_curr[,pick]
  mu1 <- c(mu.com[1] - u2*sqrt(sigma[1,1]*(1-u1)/u1),mu.com[2] - u3*sqrt(sigma[1,1]*(1-u1)/u1))
  mu2 <- c(mu.com[1] + u2*sqrt(sigma[1,1]*u1/(1-u1)),mu.com[2] + u3*sqrt(sigma[1,1]*u1/(1-u1))) 
  
  # Check if positions are within the range of the data
  if ((max(mu1[1],mu2[1]) < xup) & (min(mu1[1],mu2[1]) > xlow) & (max(mu1[2],mu2[2]) < yup) & (min(mu1[2],mu2[2]) > ylow)){
    fail_reason <- 2
    
    ##############################
    # Energy parameters
    ##############################
    # # Spectral model weight parameter
    ewt.com <- ewt[pick]
    # Approximate ewt1.range with continuously differentiable function
    ewt1.range <- c(e1_g(u1,ewt.com),e1_f(u1,ewt.com))
    #ewt1.range <- c(max(0,(ewt.com+u1-1)/u1),min(1,ewt.com/u1))
    ewt1 <- runif(1,ewt1.range[1],ewt1.range[2])
    ewt2 <- (ewt.com*w.com - ewt1*w1)/w2
    # Use the following as the spectral model Gamma weights
    pi1 <- ewt1*w1 # Called pi.11 elsewhere
    pi2 <- ewt2*w2 # Called pi.21 elsewhere 
    pi.com <- ewt.com*w.com # check = pi1 + pi2
    lambda1 <- (1-ewt1)*w1 # Called pi.12 elsewhere
    lambda2 <- (1-ewt2)*w2 # Called pi.22 elsewhere
    lambda.com <- (1-ewt.com)*w.com
    # First component
    v11 <- pi1/pi.com
    v2 <- runif(1,0,1)
    gam.com.mean1 <- eparas[pick,1]
    gam.com.shape1 <- eparas[pick,3]
    gam.mean11 <- v2*gam.com.mean1
    gam.mean21 <- (1-v11*v2)*gam.com.mean1/(1-v11)
    v4 <- rgamma(1,s,s)
    gam.shape11 <- v4*gam.com.shape1
    gam.shape21 <- pi2*gam.mean21^2/(pi.com*gam.com.mean1^2*(1+1/gam.com.shape1) - pi1*gam.mean11^2*(1+1/gam.shape11) - pi2*gam.mean21^2)
    # Second component
    v12 <- w1*(1-ewt1)/(w.com*(1-ewt.com))
    v3 <- runif(1,0,1)
    v5 <- rgamma(1,s,s)
    gam.com.mean2 <- eparas[pick,2]
    gam.com.shape2 <- eparas[pick,4]
    # We have to use gm11 as the reference to include the case gm11, gm12 < gcm1
    # Here we need to check if gm21 < gm22 because it may be that gm21 > gcm2
    gam.mean12 <- gam.mean11 + (v3/v12)*(gam.com.mean2 - gam.mean11)
    gam.mean22 <- gam.mean11 + ((1-v3)/(1-v12))*(gam.com.mean2 - gam.mean11)
    gam.shape12 <- v5*gam.com.shape2
    gam.shape22 <- lambda2*gam.mean22^2/(lambda.com*gam.com.mean2^2*(1+1/gam.com.shape2) - lambda1*gam.mean12^2*(1+1/gam.shape12) - lambda2*gam.mean22^2)
    # Checks
    #pi.com*gam.com.mean1^2*(1+1/gam.com.shape1)
    #pi1*gam.mean11^2*(1+1/gam.shape11) + pi2*gam.mean21^2*(1+1/gam.shape21)
    #lambda.com*gam.com.mean2^2*(1+1/gam.com.shape2)
    #lambda1*gam.mean12^2*(1+1/gam.shape12) + lambda2*gam.mean22^2*(1+1/gam.shape22)
    ##############################
    
    # Test of conditions
    # Check if positions are closest
    if (sum((mu1-mu2)^2) < min(apply((mu1-mu_curr[,-pick])^2,2,sum))){
      fail_reason <- 3
      # Check second small component is less than the second large component 
      # Also check that all the means are within the prior support
      if ((gam.mean21 < gam.mean22) & (sum(c(gam.mean11,gam.mean21) < emean.min)==0) & (sum(c(gam.mean12,gam.mean22) > emean.max)==0) & (min(gam.shape11,gam.shape12,gam.shape21,gam.shape22)>0)){
        fail_reason <- 4
        # Finally check if the proposed ewt1 was outside the allowed range (due to approximation)
        if ((ewt1 > max(0,(ewt.com+u1-1)/u1)) & (ewt1 < min(1,ewt.com/u1))){
          fail_reason <- NA
          
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
          A <- exp(log_accept_split_extended_full_type1(k_curr,mu_curr,mu.new,allocate_curr,new.allocate,pick,pick+1,eparas,eps.new,ewt,ewt.new,pick,w,w.new,probs.new,c(u1,u2,u3,v2,v3,v4,v5),ewt1.range))
          
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
      }
    }
  }
  
  value <- list(w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num,fail_reason,mu1,mu2,A)
  return(value)
}