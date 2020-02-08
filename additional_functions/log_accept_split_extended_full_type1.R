# Acceptance probability for split type1
log_accept_split_extended_full_type1 <- function(k,mu,newmu,allocate,newallocate,ind1,ind2,oldeps,neweps,oldew,newew,indold,wt,newwt,pupdate,us,ewt1.range){
  
  # Log complate data likelihood ratio
  like_part <- complete_log_like_extended_full(newmu,newallocate,neweps,newew,newwt,k+1)-complete_log_like_extended_full(mu,allocate,oldeps,oldew,wt,k)
  
  # Positions log prior ratio
  mupart <- log(1/img_area)
  
  # Spectral parameters log prior ratio
  if (k > 1){
    epart <- sum(log(dgamma(neweps[,3:4],ashape,arate))) - 2*log(emean.range) - sum(log(dgamma(oldeps[,3:4],ashape,arate))) 
  } else {
    epart <- sum(log(dgamma(neweps[3:4],ashape,arate))) - 2*log(emean.range) - sum(log(dgamma(oldeps[3:4],ashape,arate))) 
  }
  
  # Spectral model weight parameter log prior ratio
  ewt_part <- sum(log(dbeta(newew,ewt_ab,ewt_ab)))-sum(log(dbeta(oldew,ewt_ab,ewt_ab)))
  
  # Number log prior ratio
  kpart <- log(dpois(k+1,theta)/dpois(k,theta))
  
  # Log allocation likelihood ratio and log allocation prior probability
  count1 <- sum(newallocate[,ind1])
  count2 <- sum(newallocate[,ind2])
  wcon_part <- (wprior-1+count1)*log(newwt[ind1])+(wprior-1+count2)*log(newwt[ind2])+lgamma((k+2)*wprior)-(wprior-1+count1+count2)*log(wt[indold])-lgamma(wprior)-lgamma((k+1)*wprior)
  
  # 'Closet neighbour' rather than 'next neighbour to the right' gives these slightly more complicated move probabilities
  # (Note in both the combine and split moves the sources are ordered so source 1 has the smaller first gamma component 
  # - this meets the requirement that there is only one way to complete the proposed move)
  bk <- ifelse(k==1,1,0.25)/k # one source is chosen and split
  if (k>2){
    if (sum((newmu[,ind1]-newmu[,ind2])^2) < min(min(apply((newmu[,ind2]-mu[,-indold])^2,2,sum)),min(apply((newmu[,ind1]-mu[,-indold])^2,2,sum)))){
      dk1 <- 0.25*2/(k+1) # if mu1 and mu2 are mutually closest neighbours then there is double the probability of choosing this pair to combine
    } else {
      dk1 <- 0.25*1/(k+1)
    }
  } 
  if (k==2){
    if (sum((newmu[,ind1]-newmu[,ind2])^2) < min(sum((newmu[,ind2]-mu[,-indold])^2),sum(newmu[,ind1]-mu[,-indold])^2)){
      dk1 <- 0.25*2/(k+1) # if mu1 and mu2 are mutually closest neighbours then there is double the probability of choosing this pair to combine
    } else {
      dk1 <- 0.25*1/(k+1)
    }
  }
  if (k==1){
    dk1 <- 0.25
  }
  
  # Deal with empty sources
  indexold <- which(allocate[,indold] == 1)
  if (length(indexold)==0){
    proppart <- 0
  }
  if (length(indexold)>0){
    proppart <- -sum(log(dbinom(newallocate[indexold,ind1],1,pupdate)))
  }
  
  # Proposal density
  proppart <- log(dk1)- log(bk) + proppart - log(0.25*dbeta(us[1],u1a,u1b)*dbeta(abs(us[2]),2,2)*dbeta(abs(us[3]),2,2))
  proppart <- proppart - sum(log(dgamma(us[6:7],s,s))) # v4 and v5 
                                                       # v2 and v3 give no contribution
                                                       # t (i.e. ewt1) gives no contribution
  
  if (k > 1){
    jacobian <- log(split_jacobian((newew[ind1]-ewt1.range[1])/(ewt1.range[2]-ewt1.range[1]), oldew[indold], oldeps[indold,1], oldeps[indold,2], oldeps[indold,3], oldeps[indold,4], us[4], us[5], us[6], us[7], wt[indold], mu[1,indold],mu[2,indold], us[1], us[2], us[3]))
  } else {
    jacobian <- log(split_jacobian((newew[ind1]-ewt1.range[1])/(ewt1.range[2]-ewt1.range[1]), oldew[indold], oldeps[1], oldeps[2], oldeps[3], oldeps[4], us[4], us[5], us[6], us[7], wt[indold], mu[1],mu[2], us[1], us[2], us[3]))
  }
  
  # Parameter order: (t, ej, g.j1, g.j2, a.j1, a.j2, v2,v3, v4, v5, wj, mu.j1, mu.j2, u1, u2, u3)
  
  value <- like_part + mupart + epart + ewt_part + kpart + wcon_part + proppart + jacobian
  return(value)
}