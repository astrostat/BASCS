# King profile density PSF
# Evaluates PSF for any given center and an nx2 matrix of spatial points
psf <- function(pts,center){
  m <- length(pts)/2 
  x <- pts[1:m]-center[1]
  y <- pts[(m+1):(2*m)]-center[2]
  rad <- sqrt((x*cos(off.angle)+y*sin(off.angle))^2+(y*cos(off.angle)-x*sin(off.angle))^2/(1-ellip)^2)
  value <- psf.norm/(1+(rad/r0)^2)^slope
  return(value)
}