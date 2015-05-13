# Birth acceptance probability
log_accept_birth <- function(k,k0,news){
  
  # Positions log prior ratio
  mupart <- 0 
  
  # Number log prior ratio
  kpart <- log(dpois(k+1,theta)/dpois(k,theta))
  
  # Log allocation likelihood ratio and log allocation prior probability
  wcon_part <- (wprior-1)*log(news)+(obs_num+(k+1)*wprior-(k+1))*log(1-news)+lgamma((k+2)*wprior)-lgamma(wprior)-lgamma((k+1)*wprior)
  
  # Move probabilities
  bk <- ifelse(k==1,1,0.5)
  dk1 <- 0.5
  
  # Log proposal probability plus log Jacobian 
  propjac <- log(dk1)- log((k0+1)*bk) - log(dbeta(news,abirth,bbirth)) + (k+1)*log(1-news)
  
  value <- mupart + kpart + wcon_part + propjac
  return(value)
}