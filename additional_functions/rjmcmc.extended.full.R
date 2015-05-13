rjmcmc.extended.full <- function(mcmc_runs,online_ordering,iters,w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num){
  
  # Initialize storage
  paras <- list()
  store_log_post <- list()
  store_log_like <- list()
  
  # Variables and storage for performance checking
  count_zero <- 0
  com.prop.paras.out.of.bounds <- 0
  com1_probs <- list()
  com2_probs <- list()
  split1_probs <- list()
  split2_probs <- list()
  death_probs <- list()
  birth_probs <- list()
  birth_record <- list()
  split1_record <- list()
  split2_record <- list()
  com1_record <- list()
  com2_record <- list()
  com1_tries <- 0
  com2_tries <- 0
  split1_tries <- 0
  split2_tries <- 0
  death_tries <- 0
  birth_tries <- 0
  com1_count <- 0
  com2_count <- 0
  split1_count <- 0
  split2_count <- 0
  birth_count <- 0
  death_count <- 0
  
  for (tt in 1:iters){
    
    # Print iteration at intervals
    if (tt/print_interval == round(tt/print_interval)){
      print(paste("Iteration number: ",tt,sep=""))
      print(paste("Number of sources: ",k_curr,sep=""))
    }
    
    ##############################
    # MCMC part
    ##############################
    
    # Run `fix_runs' iterations of standard MCMC
    fixedk_update <- mcmc.extended.full(mcmc_runs,online_ordering,tt,w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num)
    new_paras <- fixedk_update[[1]]
    
    # Store parameters and log posterior (and log likelihood)
    paras[[tt]] <- new_paras
    store_log_post[[tt]] <- fixedk_update[[3]][1]
    store_log_like[[tt]] <- fixedk_update[[3]][2]
  
    # Updated parameters
    k_curr <- new_paras[1]
    mix_num <- k_curr+1
    allocate_curr <- fixedk_update[[2]]
    mu_curr <- matrix(new_paras[2:(2*k_curr+1)],2,k_curr)
    w <- new_paras[(2*(k_curr+1)):(3*k_curr+2)]
    eparas <- matrix(new_paras[(3*k_curr+3):(7*k_curr+2)],nrow=k_curr)  
    ewt <- new_paras[(7*k_curr+3):(8*k_curr+2)]
    

    ##############################
    # RJMCMC part - Update k
    ##############################
    
    ##############################
    # Split and combine move
    
    if (tt < iters){
      u <- ifelse(k_curr == 1, 1, runif(1,0,1))
      
      ##############################
      # Combine - Type 1 
      # This move is good for combining spectral distributions which are both balanced in the sense that the 
      # two Gammas components that make up the spectral distribution have similar weight. 
      ##############################
      
      if (u < 0.25){
        
        # Note for both combine moves:
        # Choose second source to be the closest one
        # Double weights sources which are 'mutually closest' - probably desirable
        # Sources are ordered so that there is only one possible way to do the move
        # Specifically, order by first gamma mean because split makes first gamma mean the smallest
        # Assumes sources are not exactly on top of each other (avoid this in initialization)
        
        split_combine_update <- combine_type1(w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num,com.prop.paras.out.of.bounds)
        
        # Record update and diagnostics
        com1_tries <- com1_tries + 1 
        com1_probs[[com1_tries]] <- split_combine_update[[11]]
        com1_record[[com1_tries]] <- cbind(split_combine_update[[9]],split_combine_update[[10]])
        com.prop.paras.out.of.bounds <- split_combine_update[[8]]
        if(split_combine_update[[6]]!=k_curr){com1_count <- com1_count + 1}
        
      }
      
      ##############################
      # Corresponding split move
      ##############################
      
      if ((u >= 0.25) & (u < 0.5)){
        
        split_combine_update <- split_type1(w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num,count_zero)
        
        # Record update and diagnostics
        split1_tries <- split1_tries + 1 
        split1_probs[[split1_tries]] <- split_combine_update[[11]]
        split1_record[[split1_tries]] <- cbind(split_combine_update[[9]],split_combine_update[[10]])
        count_zero <- split_combine_update[[10]]
        if(split_combine_update[[6]]!=k_curr){split1_count <- split1_count + 1}
        
      }
      
      ##############################
      # Combine - Type 2
      # This move is good for combining spectral distributions which are both focused 
      # on the first Gamma component
      # (this happens a lot). 
      ##############################
      
      if ((u >= 0.5) & (u < 0.75)){
        
        split_combine_update <- combine_type2(w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num,com.prop.paras.out.of.bounds)
        
        # Record update and diagnostics
        com2_tries <- com2_tries + 1 
        com2_probs[[com2_tries]] <- split_combine_update[[11]]
        com2_record[[com2_tries]] <- cbind(split_combine_update[[9]],split_combine_update[[10]])
        com.prop.paras.out.of.bounds <- split_combine_update[[8]]
        if(split_combine_update[[6]]!=k_curr){com2_count <- com2_count + 1}
        
      }
      
      ##############################
      # Corresponding split move
      ##############################
      
      if (u >= 0.75){
        
        split_combine_update <- split_type2(w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num,count_zero)
        
        # Record update and diagnostics
        split2_tries <- split2_tries + 1 
        split2_probs[[split2_tries]] <- split_combine_update[[11]]
        split2_record[[split2_tries]] <- cbind(split_combine_update[[9]],split_combine_update[[10]])
        count_zero <- split_combine_update[[10]]
        if(split_combine_update[[6]]!=k_curr){split2_count <- split2_count + 1}
        
      }
      
      
      # Update parameters
      if (split_combine_update[[6]]!=k_curr){
        w <- split_combine_update[[1]]
        allocate_curr <- split_combine_update[[2]]
        mu_curr <- split_combine_update[[3]]
        eparas <- split_combine_update[[4]]
        ewt <- split_combine_update[[5]]
        k_curr <- split_combine_update[[6]]
        mix_num <- split_combine_update[[7]]
      }
      
      #############################
      
      
      #############################
      # Birth and death move
      #############################
      
      u <- ifelse(k_curr == 1, 1, runif(1,0,1))
      
      
      #############################
      # Death
      #############################
      
      if (u < 0.5){
        
        birth_death_update <- death_move_extended_full(w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num)
        death_tries <- death_tries + 1 
        death_probs[[death_tries]] <- birth_death_update[[9]]
        if(birth_death_update[[6]]!=k_curr){death_count <- death_count + 1}
        
      }
      
      ##############################
      # Birth
      ##############################
      
      if (u >= 0.5){
        
        birth_death_update <- birth_move_extended_full(w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num)
        birth_tries <- birth_tries + 1 
        birth_record[[birth_tries]] <- birth_death_update[[8]]
        birth_probs[[birth_tries]] <- birth_death_update[[9]]
        if(birth_death_update[[6]]!=k_curr){birth_count <- birth_count + 1}
        
      }
      
      
      # Update parameters
      if (birth_death_update[[6]]!=k_curr){
        w <- birth_death_update[[1]]
        allocate_curr <- birth_death_update[[2]]
        mu_curr <- birth_death_update[[3]]
        eparas <- birth_death_update[[4]]
        ewt <- birth_death_update[[5]]
        k_curr <- birth_death_update[[6]]
        mix_num <- birth_death_update[[7]]
      }
    }
  }
  
  value <- list(paras,allocate_curr,store_log_post,store_log_like, 
                com1_probs,com2_probs,split1_probs,split2_probs,death_probs,birth_probs,
                com1_tries,com2_tries,split1_tries,split2_tries,death_tries,birth_tries,
                com1_count,com2_count,split1_count,split2_count,birth_count,death_count,birth_record,
                com1_record,com2_record,split1_record,split2_record)
  return(value)
}