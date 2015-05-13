inv.logit <- function(x){
  value <- exp(x)/(1+exp(x))
  return(value)
}