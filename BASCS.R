# BASCS algorithm: RJMCMC for separating overlapping sources
# David Jones, Harvard University

# See ``Disentangling Overlapping Astronomical Sources using Spatial and Spectral Information"
# Jones, Kashyap, van Dyk (ApJ 2015) for model and method details. 

# The allgorithm uses spatial and spectral information to find the joint posterior distribution
# of the number of sources and their spatial and spectral parameters. Specifically, the inference 
# is for the extended full model described in Jones, Kashyap, van Dyk (ApJ 2015) (this is the 
# model used in the data analysis of Section 7 of the paper and is the most detailed model we 
# considered). 

# RUN TIMES
# The algorithm is computationally intensive. 
# Example: for the small Chandra dataset (1000 photons) the full analysis takes about 30 hours 
# (200k RJMCMC iterations - which includes 200k updates of the number of sources and 2000k MCMC updates
# of the other parameters) 

home = getwd()
results <- paste(home,"/results/",sep="")

############################
# Required R packages
############################

library("MCMCpack")
library("cubature")
library("Rcpp")

#############################
# Load data
#############################

# Data - spatial and energy
# `spatial': n x 2 matrix of spatial co-ordinates of detected photons
# `energy': n x 1 vector of spectral data (PI channels)
setwd(paste(home,"/data/",sep=""))
load("small_chandra_orion_subset1522.RData")
spatial <- as.matrix(spatial)
colnames(spatial) <- NULL
rownames(spatial) <- NULL

# Basic data information
obs_num <- length(spatial[,1])
xlow <- min(spatial[,1])
xup <- max(spatial[,1])
ylow <- min(spatial[,2])
yup <- max(spatial[,2])
yl <- (yup - ylow)/2
xl <- (xup - xlow)/2
img_area <- 4*xl*yl  # Image area - rectangle assumed


###########################
# Options and settings:
###########################

# (1) Spectral model: extended_full / full / none
spectral_model <- "extended_full"

# (2) Load PSF and set parameters 
# The psf.R function has arguments `pts' (nx2 matrix of spatial points) and 
# `center' (the position at which the PSF is centered). 'psf.norm' the PSF
# normalizing constant that appears in psf.R is calculated automatically 
# below. 
setwd(paste(home,"/psf/",sep=""))
source("psf.R") # 2D King prfile density (King 1962)
psf.form <- "parametric"

if (psf.form == "parametric"){
  
  function_load_name <- ".R"
  slope <- 1.5  # power-law slope
  r0 <- 0.6  # core radius
  ellip <- 0.00573533  # ellipticity
  off.angle <- 0  # off-axis angle
  
  # Calculate normalizing constant
  source("normalize_psf.R")
  psf.norm <- 1
  psf.norm <- 1/normalize_psf(100)$integral 
  
  # Load C++ version PSF (requires psf.norm calculated above)
  sourceCpp('psf_cpp.cpp')
  
} else {

  # Alternatively load an image of the PSF 
  # User must supply image
  # List of length three:
  # psf_image[[1]] x vector
  # psf_image[[2]] y vector
  # psf_image[[3]] grid of psf values 
  function_load_name <- "_imagepsf.R"
  setwd(paste(home,"/psf/psf_image/",sep=""))
  load("psf_image.RData")
  source("psf_image_function.R") 
  source("adjust_grids.R")
  
  xpsf_grid <- psf_image[[1]] 
  ypsf_grid <- psf_image[[2]] 
  psf_image_values <- psf_image[[3]] 
  
  midpoint_alignment <- adjust_grids()
  xpsf_grid <- midpoint_alignment[[1]]
  ypsf_grid <- midpoint_alignment[[2]]
  psf_image_values <- midpoint_alignment[[3]]
  
}


# (3) Prior parameter settings
theta <- 14  # Prior: K ~ Pois(theta), where K is the numer of sources
wprior <- 1  # Prior: w ~ Dirichlet(wprior), where w is a vector giving the relative intensities of the 
# sources (and background)
ashape <- 2  # Prior: a ~ Gamma(ashape,arate), where a is the shape parameter of the spectral model
arate <- 0.5  # See above comment 
ewt_ab <- 2 # ewt ~ Beta(ewt_ab,ewt_ab), where ewt is the spectral weight parameter for the extended full model
emean.min <- min(energy)  # Prior: e ~ Uniform(emean.min,emean.max), where e is the mean parameter of the spectral model
emean.max <- max(energy)  # See above comment
emean.range <- emean.max-emean.min  # The reciporcal of the density of Uniform(emean.min,emean.max)

# (4) MCMC seetings 
adapt_end <- 10000  # Iteration number on which to stop adaptive MCMC location updates
mu_adapt_prop_sd <- 1  # Adaptive MCMC location updates (first adapt_end iterations):
# the standard deviation of the proposal distribution for 
# updating the location of source i is mu_adapt_prop_var/sqrt(ni)
# where ni is the number of photons currently assigned to 
# source i.                   
mu_fixed_prop_sd <- 0.1  # Location fixed updates (after the first adapt_end iterations):
# the standard deviation of the proposal distribution for 
# updating the location of a source is mu_fixed_prop_var
specm_sd <- 10  # Proposal: e ~ N(e_current,10*i), where e is the mean parameter of the spectral model, and
# i is the source number (in the inferred sources are ordered so that source 1 is the brightest).
# The proposal is rejected if the proposed value of e is outside (emean.min,emean.max), see (3). 
speca_sd <- 1  # Proposal: a ~ N(a_current,speca_sd), where a is the shape parameter of the spectral model
specwt_sd <- 0.1  # Proposal: ewt ~ inv.logit(rnorm(1,logit(ewt_current),specwt_sd)), where ewt is the weight parameter
# in the extended full model (whcih uses a two gamma spectral model)

# (5) RJMCMC settings
sigma <- diag(2)  # Affecte the proposal of new sources (and combination of existing sources).
# The larger the diagonal elements are the further apart two new sources are when
# an existing source is split. If the diagonal elements are very large then 
# split proposal may be very unlikely to be accepted.
u1a <- 2  # Proposal for relative intensities in split move: u1 ~ Beta(u1a,u1b), 
# w1 = u1*w, w2 = (1-u1)*w, where w1 and w2 are the relative intensities of the 
# two new sources 
u1b <- 2  # See above comment
s <- 5  # Used in proposing mean parameters of the spectral model in split proposals.
a_scale <- 20  # Same as above. The value a_scale/s should ideally be comparable to a typical shape parameter that
# we expect in the spectral model Gamma components (the smaller component is usually more important)
prop_ewta <- 10  # Proposal: ewt ~ Beta(prop_ewta,prop_ewtb), where ewt is the spectral weight parameter for the extended full 
# model. Specifically, used in type 2 split and combine move proposals. Chosen so that first component
# of the spectral model is heavily favored.
prop_ewtb <- 1 # See above. 
abirth <- 1  # Birth proposal: ewt ~ Beta(abirth,bbirth), where ewt is the relative intensity of the proposed
# new source. Our choice of the proposal parameters makes proposals very weak sources (that are 
# therefore more likely to be accepted)   
bbirth <- 400  # See above comment
# One could tune other RJMCMC parameters, but we have not done so. Neither do we claim that the values given here are 
# optimal. 

# (6) Algorithm settings
# Algorithm parameters
burnin_rjmcmc_iters <- 10^5  # Number of RJMCMC iterations for the burnin period
main_rjmcmc_iters <- 10^5 # Number of RJMCMC iterations for the main run
mcmc_runs <- 10  # Number of steps of MCMC to run between each RJMCMC proposal
print_interval <- 5000  # Number of iterations between print of current iteration number 
# and current number of sources


#############################
# Default initialization
# Users who know the approximate locations of the sources and / or 
# their relative intensities and spectral parameters can input their
# own initialization
#############################

setwd(paste(home,"/additional_functions/",sep=""))
source("initialization.R")

k_curr <- theta  # Number of sources
initialize_list <- initialization(k_curr)  

mix_num <- k_curr+1  # Number of mixture components (sources + background)
mu_curr <- initialize_list[[1]]  # Locations of sources (2 x k_curr matrix) - locations initially on edge of image
w <- initialize_list[[2]]  # Relative intensities (vector of length mix_num) 
eparas <- initialize_list[[3]]  # Shape and mean spectral parameters - disregard if spectral_model=="none"
ewt <- initialize_list[[4]]  # Spectral model weights (extended full model) - disregard unless spectral_model=="extended_full"
allocate_curr <- initialize_list[[5]]  # Allocation of photons to sources (and background) - initially obs_num x mix_num matrix of zeros


#############################
# Code required
#############################

setwd(paste(home,"/additional_functions/",sep=""))
# MCMC 
source(paste("mcmc.extended.full",function_load_name,sep=""))
source("spectral_post.R")
source(paste("extended_full_log_posterior",function_load_name,sep=""))
source("logit.R")
source("inv.logit.R")

# RJMCMC
source(paste("combine_type1",function_load_name,sep=""))
source(paste("split_type1",function_load_name,sep=""))
source(paste("combine_type2",function_load_name,sep=""))
source(paste("split_type2",function_load_name,sep=""))
source("birth_move_extended_full.R")
source("death_move_extended_full.R")
source(paste("complete_log_like_extended_full",function_load_name,sep=""))
source("log_accept_split_extended_full_type1.R")
source("log_accept_split_extended_full_type2.R")
source("rjmcmc.extended.full.R")
source("split_jacobian.R")
source("log_accept_birth.R")
source('e1_f.R')
source('e1_g.R')


#############################
# Initial RJMCMC run
#############################

online_ordering <- "none"  # How the sources should be ordered after each iteration. Options: "none" / "reference"
# The option "reference" means that the sources will be ordered by their suspected
# relative intensities. In particular, the current sources will be matched 
# (as near as possible) to reference positions. The reference positions are typically
# obtained by running the RJMCMC with no ordering until convergence and then using 
# the final iteration as the reference for new draws. A reference works well if the source
# position posteriors do not have substantial overlap. 

# Record start time
start <- proc.time()

initial_run <- rjmcmc.extended.full(mcmc_runs,online_ordering,burnin_rjmcmc_iters,w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num)

# Extract final parameters and allocation matrix
new_paras <- initial_run[[1]][[burnin_rjmcmc_iters]]
k_curr <- new_paras[1]
mix_num <- k_curr+1
mu_curr <- matrix(new_paras[2:(2*k_curr+1)],2,k_curr)
w <- new_paras[(2*(k_curr+1)):(3*k_curr+2)]
eparas <- matrix(new_paras[(3*k_curr+3):(7*k_curr+2)],nrow=k_curr)  
ewt <- new_paras[(7*k_curr+3):(8*k_curr+2)]
allocate_curr <- initial_run[[2]]


#############################
# Main RJMCMC run
#############################

online_ordering <- "reference"
no_guess <- k_curr
mu_guess <- mu_curr[,order(w[-mix_num],decreasing=TRUE)]
adapt_end <- 0  # Iteration number on which to stop adaptive MCMC location updates

main_run <- rjmcmc.extended.full(mcmc_runs,online_ordering,main_rjmcmc_iters,w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num)

# Extract final parameters and allocation matrix
new_paras <- main_run[[1]][[main_rjmcmc_iters]]
k_curr <- new_paras[1]
mix_num <- k_curr+1
mu_curr <- matrix(new_paras[2:(2*k_curr+1)],2,k_curr)
w <- new_paras[(2*(k_curr+1)):(3*k_curr+2)]
eparas <- matrix(new_paras[(3*k_curr+3):(7*k_curr+2)],nrow=k_curr)  
ewt <- new_paras[(7*k_curr+3):(8*k_curr+2)]
allocate_curr <- main_run[[2]]

# Save output
save(list=ls(),file=paste(results,"small_chandra_bascs_analysis.RData",sep=""))

# Show time taken
end <- proc.time()
print(end-start)

# Simple plot of final positions
plot(spatial)
points(t(mu_curr),pch=16,col=2)


