# Acceptance probability for split type2
log_accept_split_extended_full_type2 <- function(k,mu,newmu,allocate,newallocate,ind1,ind2,oldeps,neweps,oldew,newew,indold,wt,newwt,pupdate,us){
  
  # Log complate data likelihood ratio
  like_part <- complete_log_like_extended_full(newmu,newallocate,neweps,newew,newwt,k+1)-complete_log_like_extended_full(mu,allocate,oldeps,oldew,wt,k)
  
  # Positions log prior ratio
  mupart <-log(1/img_area)
  
  # Spectral parameters log prior ratio
  epart <- sum(log(dgamma(neweps[,3:4],ashape,arate))) - 2*log(emean.range) - sum(log(dgamma(oldeps[,3:4],ashape,arate))) # means are uniformly distributed across the range of the data
  
  # Spectral model weight parameter log prior ratio
  ewt_part <- sum(log(dbeta(newew,ewt_ab,ewt_ab)))-sum(log(dbeta(oldew,ewt_ab,ewt_ab)))
  
  # Number log prior ratio
  kpart <- log(dpois(k+1,theta)/dpois(k,theta))
  
  # Log allocation likelihood ratio and log allocation prior probability
  count1 <- sum(newallocate[,ind1])
  count2 <- sum(newallocate[,ind2])
  wcon_part <- (wprior-1+count1)*log(newwt[ind1])+(wprior-1+count2)*log(newwt[ind2])+lgamma((k+2)*wprior)-(wprior-1+count1+count2)*log(wt[indold])-lgamma(wprior)-lgamma((k+1)*wprior)
  
  # 'Closet neighbour' rather than 'next neighbour to the right' gives these slightly more complicated move probabilities
  # (Note in both the combine and eplit moves the sources are ordered so source 1 has the smaller first gamma component 
  # - this meets the requirement that there is only one way to complete the proposed move)
  bk <- ifelse(k==1,1,0.25)/k # one source is chosen and split
  if (sum((newmu[,ind1]-newmu[,ind2])^2) < min(min(apply((newmu[,ind2]-mu[,-indold])^2,2,sum)),min(apply((newmu[,ind1]-mu[,-indold])^2,2,sum)))){
    dk1 <- 0.25*2/(k+1) # if mu1 and mu2 are mutually closest neighbours then there is double the probability of choosing this pair to combine
  } else {
    dk1 <- 0.25*1/(k+1)
  }
  
  # Deal with empty sources
  proppart <- 0
  indexold <- which(allocate[,indold] == 1)
  if (length(indexold)>0){
    proppart <- -sum(log(dbinom(newallocate[indexold,ind1],1,pupdate)))
  }
  
  # Proposal density
  proppart <- log(dk1)- log(bk) + proppart - log(0.25*dbeta(abs(us[2]),2,2)*dbeta(abs(us[3]),2,2))
  proppart <- proppart - log(dbeta(us[4],prop_ewta,prop_ewtb)) - log(dbeta(us[5],prop_ewta,prop_ewtb))  # ewt1 and ewt2
  proppart <- proppart - log(dbeta(us[6],1,s)) - log(dbeta(us[7],1,s))  # alpha12/20 and alpha22/20
  
  # Jacobian
  jacobian <- log(abs(a_scale^2*(emean.max-mu[1,indold])*(emean.max-mu[2,indold])*wt[indold]*sigma[1,1]/(us[1]*(1-us[1]))))
  
  value <- like_part + mupart + epart + ewt_part + kpart + wcon_part + proppart + jacobian
  return(value)
}