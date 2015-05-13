# Extended full model death move
death_move_extended_full <- function(w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num){
  
  # Can only kill empty sources
  A <- 0
  pick <- 0
  count_short <- apply(allocate_curr[,1:k_curr],2,sum)
  if (min(count_short) == 0){
    index_short <- which(count_short==0)
    pick <- ifelse(length(index_short) > 1,sample(index_short,1),index_short)
    
    # Acceptance probability
    count_short_new <- apply(allocate_curr[,1:k_curr],2,sum)
    empty <- sum(count_short_new == 0)
    A <- exp(-log_accept_birth(k_curr-1,empty-1,w[pick]))
    
    utest <- runif(1,0,1)
    if (is.na(A)==0){
      if (utest < A){
        mu_curr <- matrix(mu_curr[,-pick],nrow=2,ncol=k_curr-1)
        eparas <- matrix(eparas[-pick,],nrow=k_curr-1,ncol=4)
        ewt <- ewt[-pick]
        k_curr <- k_curr-1
        mix_num <- k_curr+1
        w <- w[-pick]/sum(w[-pick])
        allocate_curr <- allocate_curr[,-pick]
      }
    }
  }
  
  value <- list(w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num,pick,A)
  return(value)
}