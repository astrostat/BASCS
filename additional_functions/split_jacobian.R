split_jacobian <- function(t, ej, g.j1, g.j2, a.j1, a.j2, v2,v3, v4, v5, wj, mu.j1, mu.j2, u1, u2, u3){

#####################
# Initial matrix
#####################

jab.mat <- matrix(NA,10,16)

colnames(jab.mat) <- c("t","ej","g.j1","g.j2","a.j1","a.j2","v2","v3","v4","v5","wj","mu.j1","mu.j2","u1","u2","u3")

rownames(jab.mat) <- c("ej1","ej2","g.j11","g.j21","g.j12","g.j22","a.j11","a.j21","a.j12","a.j22") #,"w1","w2","mu.j11","mu.j12","mu.j21","mu.j22")


###################################################################
# Eparas part
###################################################################


#####################
# Define variables
#####################

ej1 <- t*e1_f(u1,ej) + (1-t)*e1_g(u1,ej)

wj1 <- u1*wj

wj2 <- (1-u1)*wj

pi.j11 <- ej1*wj1

pi.j1 <- ej*wj

v11 <- pi.j11/pi.j1

pi.j21 <- (1-v11)*pi.j1

v12 <- u1*(1-ej1)/(1-ej)

pi.j2 <- (1-ej)*wj

pi.j12 <- v12*pi.j2

pi.j22 <- (1-v12)*pi.j2

ej2 <- (ej-u1*ej1)/(1-u1)  # Check: pi.j2/wj2

g.j11 <- v2*g.j1
  
g.j21 <- (1-v11*v2)/(1-v11)*g.j1

g.j12 <- g.j11 + (v3/v11)*(g.j2 - g.j11)

g.j22 <- g.j11 + ((1-v3)/(1-v11))*(g.j2 - g.j11)

a.j11 <- v4*a.j1

a.j12 <- v5*a.j2


#####################
# e1.f(u_1,ej)
#####################

if (ej/u1 < 1){
  df_du1 <- -ej/u1^2
} else {
  df_du1 <- (10/ej)*exp(10*(u1/ej-1))*log(ej/u1) - (1/u1)*exp(10*(u1/ej-1))
}

if (ej/u1 < 1){
  df_dej <- 1/u1
} else {
  df_dej <- (1/ej)*exp(10*(u1/ej-1)) - (10*u1/ej^2)*exp(10*(u1/ej-1))*log(ej/u1)
}


######################
# e1.g(u_1,ej)
######################

if ((ej+u1-1)/u1 > 0){
  dg_du1 <- (1-ej)/u1^2
} else {
  dg_du1 <- ((1-ej)/u1^2)*exp(10*((ej+u1-1)/u1)) + 10*((ej+u1-1)/u1)*((1-ej)/u1^2)*exp(10*((ej+u1-1)/u1))
}

if ((ej+u1-1)/u1 > 0){
  dg_dej <- 1/u1
} else {
  dg_dej <- (1/u1)*exp(10*((ej+u1-1)/u1)) + (10/u1)*((ej+u1-1)/u1)*exp(10*((ej+u1-1)/u1))
}


######################
# ej1
######################

dej1_dt <- e1_f(u1,ej)-e1_g(u1,ej)

dej1_du1 <- t*df_du1 + (1-t)*dg_du1

dej1_dej <- t*df_dej + (1-t)*dg_dej

#colnames(jab.mat)[c(1,2,14)]
#rownames(jab.mat)[1]
jab.mat[1,c(1,2,14)] <- c(dej1_dt, dej1_dej, dej1_du1)

jab.mat[1,-c(1,2,14)] <- 0


######################
# ej2
######################

dej2_dt <- -(u1*dej1_dt)/(1-u1)

dej2_du1 <- (ej-u1*ej1)/(1-u1)^2 - (ej1+u1*dej1_du1)/(1-u1)

dej2_dej <- (1-u1*dej1_dej)/(1-u1)

#rownames(jab.mat)[2]
jab.mat[2,c(1,2,14)] <- c(dej2_dt, dej2_dej, dej2_du1)

jab.mat[2,-c(1,2,14)] <- 0


######################
# v11 = u1*ej1/ej
######################

dv11_du1 <- ej1/ej + (u1*dej1_du1)/ej

dv11_dej <- -(u1*ej1)/ej^2 + (u1*dej1_dej)/ej
 
dv11_dt <- (u1*dej1_dt)/ej


######################
# v12 = u1*(1-ej1)/(1-ej)
######################

dv12_du1 <- (1-ej1)/(1-ej) - (u1*dej1_du1)/(1-ej)

dv12_dej <- (u1*(1-ej1))/(1-ej)^2 - (u1*dej1_dej)/(1-ej)

dv12_dt <- -(u1*dej1_dt)/(1-ej)


######################
# g.j11
######################

dg.j11_dv2 <- g.j1

dg.j11_dg.j1 <- v2

#colnames(jab.mat)[c(3,7)]
#rownames(jab.mat)[3]
jab.mat[3,c(3,7)] <- c(dg.j11_dg.j1, dg.j11_dv2)

jab.mat[3,-c(3,7)] <- 0


######################
# g.j21
######################

dg.j21_term <- (-(v2/(1-v11)) + ((1-v11*v2)/(1-v11)^2))*g.j1

dg.j21_dt <- dv11_dt*dg.j21_term

dg.j21_du1 <- dv11_du1*dg.j21_term

dg.j21_dej <- dv11_dej*dg.j21_term

dg.j21_dv2 <- -(v11/(1-v11))*g.j1

dg.j21_dg.j1 <- (1-v11*v2)/(1-v11)

#colnames(jab.mat)[c(1,2,3,7,14)]
#rownames(jab.mat)[4]
jab.mat[4,c(1,2,3,7,14)] <- c(dg.j21_dt, dg.j21_dej, dg.j21_dg.j1, dg.j21_dv2, dg.j21_du1)

jab.mat[4,-c(1,2,3,7,14)] <- 0


######################
# g.j12
######################

dg.j12_dg.j2 <- v3/v12

dg.j12_dg.j1 <- v2 - (v3/v12)*v2

dg.j12_dv2 <- g.j1 - (v3/v12)*g.j1

dg.j12_dv3 <- (1/v12)*(g.j2 - g.j11)

dg.j12_term <- (-(v3/v12^2)*(g.j2 - g.j11))
  
dg.j12_dt <- dv12_dt*dg.j12_term

dg.j12_du1 <- dv12_du1*dg.j12_term

dg.j12_dej <- dv12_dej*dg.j12_term

#colnames(jab.mat)[c(1,2,3,4,7,8,14)]
#rownames(jab.mat)[5]
jab.mat[5,c(1,2,3,4,7,8,14)] <- c(dg.j12_dt, dg.j12_dej, dg.j12_dg.j1, dg.j12_dg.j2, dg.j12_dv2, dg.j12_dv3, dg.j12_du1)

jab.mat[5,-c(1,2,3,4,7,8,14)] <- 0


######################
# g.j22
######################

dg.j22_dg.j1 <- v2 + ((1-v3)/(1-v12))*(-v2)

dg.j22_dg.j2 <- (1-v3)/(1-v12)

dg.j22_dv2 <- g.j1 + ((1-v3)/(1-v12))*(-g.j1)

dg.j22_dv3 <- - (1/(1-v12))*(g.j2 - g.j11)

dg.j22_term <- ((1-v3)/(1-v12)^2)*(g.j2 - g.j11)

dg.j22_dt <- dv12_dt*dg.j22_term

dg.j22_du1 <- dv12_du1*dg.j22_term

dg.j22_dej <- dv12_dej*dg.j22_term

#colnames(jab.mat)[c(1,2,3,4,7,8,14)]
#rownames(jab.mat)[6]
jab.mat[6,c(1,2,3,4,7,8,14)] <- c(dg.j22_dt, dg.j22_dej, dg.j22_dg.j1, dg.j22_dg.j2, dg.j22_dv2, dg.j22_dv3, dg.j22_du1)

jab.mat[6,-c(1,2,3,4,7,8,14)] <- 0


######################
# a.j11
######################

da.j11_dv4 <- a.j1

da.j11_da.j1 <- v4

#colnames(jab.mat)[c(5,9)]
#rownames(jab.mat)[7]
jab.mat[7,c(5,9)] <- c(da.j11_da.j1, da.j11_dv4)

jab.mat[7,-c(5,9)] <- 0


######################
# a.j12
######################

da.j12_dv5 <- a.j2

da.j12_da.j2 <- v5

#colnames(jab.mat)[c(6,10)]
#rownames(jab.mat)[9]
jab.mat[9,c(6,10)] <- c(da.j12_da.j2, da.j12_dv5)

jab.mat[9,-c(6,10)] <- 0


######################
# pi.j11
######################

dpi.j11_dwj <- v11*ej

dpi.j11_dej <- v11*wj + dv11_dej*wj*ej

dpi.j11_du1 <- wj*ej*dv11_du1

dpi.j11_dt <- wj*ej*dv11_dt


######################
# pi.j21
######################

dpi.j21_dwj <- (1-v11)*ej

dpi.j21_dej <- (1-v11)*wj - dv11_dej*wj*ej

dpi.j21_du1 <- -wj*ej*dv11_du1

dpi.j21_dt <- -wj*ej*dv11_dt


######################
# pi.j12
######################

dpi.j12_dwj <- v12*ej

dpi.j12_dej <- v12*wj + dv12_dej*wj*ej

dpi.j12_du1 <- wj*ej*dv12_du1

dpi.j12_dt <- wj*ej*dv12_dt


######################
# pi.j22
######################

dpi.j22_dwj <- (1-v12)*ej

dpi.j22_dej <- (1-v12)*wj - dv12_dej*wj*ej

dpi.j22_du1 <- -wj*ej*dv12_du1

dpi.j22_dt <- -wj*ej*dv12_dt


######################
# a.j21
######################

# Variables other than g.j1, a.j1, ej, wj which matter are:
# t, u1, v2, v4

A.1 <- pi.j21*g.j21^2

B.1 <- pi.j1*g.j1^2*(1+1/a.j1) - pi.j11*g.j11^2*(1+1/a.j11) - pi.j21*g.j21^2
  
derivs1 <- matrix(NA,8,5)

# g.j21, g.j11, a.j11
# t, v2, v4, u1, wj, ej, g.j1, a.j1
#colnames(jab.mat)[c(1,7,9,14,11,2,3,5)]
#rownames(jab.mat)[c(4,3,7)]
derivs1[,c(2,4,5)] <- t(jab.mat[c(4,3,7),c(1,7,9,14,11,2,3,5)])

# pi.j11
derivs1[,3] <- c(dpi.j11_dt, 0, 0, dpi.j11_du1, dpi.j11_dwj, dpi.j11_dej, 0, 0)

# pi.j21
derivs1[,1] <- c(dpi.j21_dt, 0, 0, dpi.j21_du1, dpi.j21_dwj, dpi.j21_dej, 0, 0)

da.j21_derivs <- matrix(NA,8,1)

for (di in 1:8){
  
  term1 <- (1/B.1)*(derivs1[di,1]*g.j21^2 + 2*pi.j21*derivs1[di,2]*g.j21) 

  term2 <- -(A.1/B.1^2)*(-derivs1[di,3]*g.j11^2*(1+1/a.j11) - 2*pi.j11*derivs1[di,4]*g.j11*(1+1/a.j11)) 

  term3 <- -(A.1/B.1^2)*(-pi.j11*g.j11^2*(-1/a.j11^2*derivs1[di,5])) 

  term4 <- -(A.1/B.1^2)*(-derivs1[di,1]*g.j21^2 - 2*pi.j21*derivs1[di,2]*g.j21) 

  da.j21_derivs[di] <- term1 + term2 + term3 + term4
  
}

#colnames(jab.mat)[c(1,7,9,14)]
#rownames(jab.mat)[8]
jab.mat[8,c(1,7,9,14)] <- da.j21_derivs[1:4]

#colnames(jab.mat)[c(11,2,3,5)]
da.j21_dwj <- da.j21_derivs[5] -(A.1/B.1^2)*(ej*g.j1^2*(1+1/a.j1)) # This is zero because the terms cancel

da.j21_dej <- da.j21_derivs[6] -(A.1/B.1^2)*(wj*g.j1^2*(1+1/a.j1)) 

da.j21_dg.j1 <- da.j21_derivs[7] -(A.1/B.1^2)*(2*pi.j1*g.j1*(1+1/a.j1)) 

da.j21_da.j1 <- da.j21_derivs[8] -(A.1/B.1^2)*(pi.j1*g.j1^2*(-1/a.j1^2)) 

da.j21_derivs[5:8] <- c(da.j21_dwj, da.j21_dej, da.j21_dg.j1 ,da.j21_da.j1)

jab.mat[8,c(11,2,3,5)] <- da.j21_derivs[5:8]

jab.mat[8,-c(1,7,9,14,11,2,3,5)] <- 0


######################
# a.j22
######################

# Variables other than g.j2, a.j2, ej, wj which matter are:
# t, u1, v3, v5

A.2 <- pi.j22*g.j22^2

B.2 <- pi.j2*g.j2^2*(1+1/a.j2) - pi.j12*g.j12^2*(1+1/a.j12) - pi.j22*g.j22^2

derivs2 <- matrix(NA,8,5)

#      g.j22,     g.j12, a.j12
# t, v3, v5, u1, wj, ej, g.j2, a.j2
#colnames(jab.mat)[c(1,8,10,14,11,2,4,6)]
#rownames(jab.mat)[c(6,5,9)]
derivs2[,c(2,4,5)] <- t(jab.mat[c(6,5,9),c(1,8,10,14,11,2,4,6)])

# pi.j11
derivs2[,3] <- c(dpi.j12_dt, 0, 0, dpi.j12_du1, dpi.j12_dwj, dpi.j12_dej, 0, 0)

# pi.j21
derivs2[,1] <- c(dpi.j22_dt, 0, 0, dpi.j22_du1, dpi.j22_dwj, dpi.j22_dej, 0, 0)

da.j22_derivs <- matrix(NA,8,1)

for (di in 1:8){
  
  term1 <- (1/B.2)*(derivs2[di,1]*g.j22^2 + 2*pi.j22*derivs2[di,2]*g.j22) 
  
  term2 <- -(A.2/B.2^2)*(-derivs2[di,3]*g.j12^2*(1+1/a.j12) - 2*pi.j12*derivs2[di,4]*g.j12*(1+1/a.j12)) 
  
  term3 <- -(A.2/B.2^2)*(-pi.j12*g.j12^2*(-1/a.j12^2*derivs2[di,5])) 
  
  term4 <- -(A.2/B.2^2)*(-derivs2[di,1]*g.j22^2 - 2*pi.j22*derivs2[di,2]*g.j22) 
  
  da.j22_derivs[di] <- term1 + term2 + term3 + term4
  
}

#colnames(jab.mat)[c(1,8,10,14)]
#rownames(jab.mat)[10]
jab.mat[10,c(1,8,10,14)] <- da.j22_derivs[1:4]

da.j22_dwj <- da.j22_derivs[5] -(A.2/B.2^2)*((1-ej)*g.j2^2*(1+1/a.j2)) 

da.j22_dej <- da.j22_derivs[6] -(A.2/B.2^2)*(-wj*g.j2^2*(1+1/a.j2)) 

da.j22_dg.j2 <- da.j22_derivs[7] -(A.2/B.2^2)*(2*pi.j2*g.j2*(1+1/a.j2)) 

da.j22_da.j2 <- da.j22_derivs[8] -(A.2/B.2^2)*(pi.j2*g.j2^2*(-1/a.j2^2)) 

da.j22_derivs[5:8] <- c(da.j22_dwj, da.j22_dej, da.j22_dg.j2 ,da.j22_da.j2)

#colnames(jab.mat)[c(11,2,4,6)]
jab.mat[10,c(11,2,4,6)] <- da.j22_derivs[5:8]

jab.mat[10,-c(1,8,10,14,11,2,4,6)] <- 0


######################
# Compute determinant for full transformation including spatial position part
######################

value <- abs(det(jab.mat[1:10,1:10])*wj*sigma[1,1]/(u1*(1-u1))) # The covariance matrix sigma gives the VARIANCE part

return(value)

}
