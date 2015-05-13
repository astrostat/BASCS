normalize_psf <- function(bound){

  norm_fun <- function(pts){
    return(psf(pts,c(0,0)))
  }
  
  value <- adaptIntegrate(norm_fun,c(-bound,-bound),c(bound,bound))
  return(value)
}