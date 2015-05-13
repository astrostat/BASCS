e1_f <- function(u1.var,ej.var){
  if (ej.var/u1.var < 1){
    value <- ej.var/u1.var
  } else {
    value <- 1 + exp(10*(u1.var/ej.var-1))*log(ej.var/u1.var)
  }
  return(value)
}