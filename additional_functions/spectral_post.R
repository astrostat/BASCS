spectral_post <- function(m1,a1,m2,a2,ew,dat){
  # Proposals are all symmetric and so can be ignored
  priorm <- -2*log(emean.range)
  priora <- log(dgamma(a1,ashape,arate))+log(dgamma(a2,ashape,arate))
  prior_ew <- log(dbeta(ew,ewt_ab,ewt_ab))
  value <- sum(log(ew*dgamma(dat,a1,a1/m1)+(1-ew)*dgamma(dat,a2,a2/m2)))+priorm+priora+prior_ew
  return(value)
}