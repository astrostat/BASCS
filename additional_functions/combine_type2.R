# Extended full model combine move type 2 
combine_type2 <- function(w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num,com.prop.paras.out.of.bounds){ 
  
  # Initializae acceptance probability
  A <- 0
  
  # Pick sources to combine and order them by the their spectral component with lower mean
  # (needed to ensure uniqueness of move)
  pick.unordered <- sample(1:k_curr,1)
  mu1.unordered <- mu_curr[,pick.unordered]
  pick2.unordered <- order(apply((mu_curr-mu1.unordered)^2,2,sum))[2]
  if (order(c(eparas[pick.unordered,1],eparas[pick2.unordered,1]))[1] == 1){
    pick <- pick.unordered
    pick2 <- pick2.unordered
  } else {
    pick <- pick2.unordered
    pick2 <- pick.unordered
  }
  mu1 <- mu_curr[,pick]
  mu2 <- mu_curr[,pick2]
  w1 <- w[pick]
  w2 <- w[pick2]
  w.com <- w1+w2
  mu.com <- (w1*mu1+w2*mu2)/w.com
  
  ##############################
  # Energy parameters 
  ewt1 <- ewt[pick]
  ewt2 <- ewt[pick2]
  # In ideal situations for move type 2, ewt.com should be about 0.5 and ewt1, ewt2 about 1
  ewt.com <- w1/w.com 
  # Parameters
  gam.mean11 <- eparas[pick,1]
  gam.mean21 <- eparas[pick2,1]
  gam.shape11 <- eparas[pick,3]
  gam.shape21 <- eparas[pick2,3]
  gam.mean12 <- eparas[pick,2]
  gam.mean22 <- eparas[pick2,2]
  gam.shape12 <- eparas[pick,4]
  gam.shape22 <- eparas[pick2,4]
  # Combine: here we simply ignore the second components of the spectral distributions as we are trying to 
  # deal with the case where the first components have nearly all the weight. 
  # First component
  gam.com.mean1 <- gam.mean11
  gam.com.shape1 <- gam.shape11
  # Second component
  gam.com.mean2 <- gam.mean21
  gam.com.shape2 <- gam.shape21
  ##############################
  
  # Define proposal parameter vectors
  # Positions
  mu.interm <- mu_curr
  mu.interm[,pick] <- mu.com
  mu.new <- mu.interm[,-pick2]
  # Weights
  w.interm <- w
  w.interm[pick] <- w.com
  w.new <- w.interm[-pick2]
  # Free parameters for reverse split move
  u1 <- w1/(w1+w2)
  u2 <- (mu.com[1]-mu1[1])/(sqrt(sigma[1,1]*(1-u1)/u1))
  u3 <- (mu.com[2]-mu1[2])/(sqrt(sigma[1,1]*(1-u1)/u1))
  # Energy paras
  ewt.interm <- ewt
  ewt.interm[pick] <- ewt.com
  ewt.new <- ewt.interm[-pick2]
  eps.interm <- eparas
  eps.interm[pick,] <- c(gam.com.mean1,gam.com.mean2,gam.com.shape1,gam.com.shape2)
  eps.new <- eps.interm[-pick2,]
  # Free parameters for reverse split move
  v2 <- (gam.mean12-gam.com.mean1)/(emean.max-gam.com.mean1)
  v3 <- (gam.mean22-gam.com.mean2)/(emean.max-gam.com.mean2)
  v4 <- gam.shape12/a_scale
  v5 <- gam.shape22/a_scale
  
  # Check u2, u3 are between 0 and 1
  if ((max(abs(u2),abs(u3)) < 1)){
    # Allocation
    allocate.interm <- allocate_curr
    allocate.interm[,pick] <- allocate_curr[,pick]+allocate_curr[,pick2]
    allocate.new <- allocate.interm[,-pick2]
    # Check 
    # apply(allocate.new,2,sum)
    indcom <- which((allocate_curr[,pick]+allocate_curr[,pick2])==1)
    p1 <- w1*psf_cpp(off.angle,ellip,slope,psf.norm,r0,spatial[indcom,],mu1)*(ewt1*dgamma(energy[indcom],gam.shape11,gam.shape11/gam.mean11)+(1-ewt1)*dgamma(energy[indcom],gam.shape12,gam.shape12/gam.mean12))
    p2 <- w2*psf_cpp(off.angle,ellip,slope,psf.norm,r0,spatial[indcom,],mu2)*(ewt2*dgamma(energy[indcom],gam.shape21,gam.shape21/gam.mean21)+(1-ewt2)*dgamma(energy[indcom],gam.shape22,gam.shape22/gam.mean22))
    probs.old <- p1/(p1+p2)
    # sum(allocate_curr[indcom,pick]*log(probs.old)) + sum(allocate_curr[indcom,pick2]*log(1-probs.old))
    # sum(allocate_curr[indcom,pick]*log(probs[pick,indcom])) + sum(allocate_curr[indcom,pick2]*log(1-probs[pick2,indcom]))
    # Acceptance probability
    split.index <- which(eps.new[,1]==gam.com.mean1)
    A <- exp(-log_accept_split_extended_full_type2(k_curr-1,mu.new,mu_curr,allocate.new,allocate_curr,pick,pick2,eps.new,eparas,ewt.new,ewt,split.index,w.new,w,probs.old,c(ewt.com,u2,u3,ewt1,ewt2,v4,v5)))
    
    utest <- runif(1,0,1)
    if (is.na(A)==0){
      if (utest < A){
        w <- w.new
        allocate_curr <- allocate.new
        mu_curr <- mu.new
        eparas <- eps.new
        ewt <- ewt.new
        k_curr <- k_curr-1
        mix_num <- k_curr +1
      }
    }
  } else {
    com.prop.paras.out.of.bounds <- com.prop.paras.out.of.bounds + 1
  }
  
  value <- list(w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num,com.prop.paras.out.of.bounds,mu1,mu2,A)
  return(value)
}