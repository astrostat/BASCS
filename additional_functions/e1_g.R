e1_g <- function(u1.var,ej.var){
  if ((ej.var+u1.var-1)/u1.var > 0){
    value <- (ej.var+u1.var-1)/u1.var
  } else {
    value <- ((ej.var+u1.var-1)/u1.var)*exp(10*((ej.var+u1.var-1)/u1.var))
  }
  return(value)
}