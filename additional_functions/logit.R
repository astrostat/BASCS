logit <- function(x){
  value <- log(x/(1-x))
  return(value)
}