#####################################
# Authors: Alex Tarr
# Created: 2015/01/09
# Last Revision: N/A
# Institution: Princeton University
# Purpose: Optimize temperatures for 
#          parallel tempering
#####################################

## Inputs
tBeta <- 400 ##Edit this line for target beta of state
betaMin <- 0.01 ##Smallest beta assumed to give good mixing

## Set params (for running swMH)
params <- expand.grid(state = "ms",
                      eprob = 0.05, marginpct = 1,
                      lambda = 18, pnum = 1,
                      initbeta = 0,
                      initbetadiss = 0,
                      initbetaswitch = 0,
                      annealbeta = 0, annealbetadiss = 0,
                      annealbetapop = 1, annealbetaswitch = 0,
                      targbetapop = 0,
                      bybetapop = 0,
                      weightpow = 0,
                      thin = 1,
                      wd = "/scratch/network/bfifield/segregation/data/",
                      logdir = "/scratch/network/bfifield/segregation/code/slurm/",
                      dwd = "/scratch/network/bfifield/segregation/data/simRuns/",
                      codedir = "/scratch/network/bfifield/segregation/code/redist-pkg/")

## Modify pnum
params$pnum <- 1:2

## Set working directory | state
setwd(paste(params$wd[1], params$state[1], sep = ""))

## Load data
load(paste(getwd(), "/algdat.RData", sep = ""))

###################
# Initializations #
###################

#Convergence flag
converge <- FALSE

#Initial districting plan
if(substr(state, 1, 7) == "testset"){
  cdsAcc <- cdmat[,sample(1:ncol(cdmat), 1)]
  cdsProp <- cdmat[,sample(1:ncol(cdmat), 1)]
} else{
  cdsAcc <- eval(parse(text = paste("geodat$cds", 1, sep = "")))
  cdsProp <- eval(parse(text = paste("geodat$cds", 2, sep = "")))
}

#Scaling for determining proposed beta
initRho <- -1.5
rho <- initRho
#Chain temperatures
betaseq <- -tBeta*c(1,1/(1+exp(rho)))
#Temperature index
i <- 1
#Iteration index
n <- 1

## Function for generating samples
ecuts <- function(cds,params,betapop,al.pc,geodat){
  ## Load packages
  library("redist"); library("BARD"); library("maptools"); library("methods")
  
  ##########################
  ## Parameter Extraction ##
  ##########################
  
  ## Set State
  state <- params$state
  
  ## Set directories
  dwd <- params$dwd
  codedir <- params$codedir
  logdir <- params$logdir
  
  ## Determine district type
  if(substr(state, 1, 7) == "testset"){
    dists <- length(unique(cdmat[,1]))
  } 
  else{
    dists <- length(unique(geodat$cds))
  }
  
  ## Edgecut prob. (prob. of turning edge off = 1 - eprob)
  eprob <- params$eprob
  
  ## Population constraint
  margin.pct <- params$marginpct
  
  ## Lambda
  lambda <- params$lambda
  
  ## Par num
  pnum <- params$pnum
  
  ## Compactness beta
  beta <- params$initbeta
  betadiss <- params$initbetadiss
  betaswitch <- params$initbetaswitch

  ## Thinning amount
  thin <-  params$thin
  
  ## Weights power
  wpow <- params$weightpow
  
  ## Pop. constraints quantities
  parity <- sum(geodat$pop) / dists
  margin <- round(parity * margin.pct)
  
  ## Empty beta vector and weights for simulated tempering
  bvec <- rep(1, 5)  
  betaweights <- rep(1, length(bvec))
  
  #########################
  ## Run the simulations ##
  #########################
  
  samp <- swMH(al.pc, cds, cds, 1, eprob,
               geodat$pop, geodat$blackhisp, parity, margin,
               dists, lambda, ssdmat,
               beta = beta, betadiss = betadiss, betapop = betapop,
               betaswitch = betaswitch,
               betavec = bvec, betadissvec = bvec,
               betapopvec = bvec, betaswitchvec = bvec,
               betaweights = betaweights)
  
  return(samp)
}

## Iterative algorithm for choosing betaseq
while(!converge){
  
  ## Generate samples
  #Accepted temperature
  ecutsAcc <- ecuts(cdsAcc,params,betaseq[i],al.pc,geodat)
  #Proposed adjacent temperature
  ecutsProp <- ecuts(cdsProp,params,betaseq[i+1],al.pc,geodat)
  
  ## Update current district
  cdsAcc <- ecutsAcc[[1]][,1]
  cdsProp <- ecutsProp[[1]][,1]
  
  ## Get likelihoods
  likePop.Acc <- ecutsAcc[[18]]
  likePop.Prop <- ecutsProp[[18]]
  
  ## Compute acceptance probability
  alpha <- min(1,exp((betaseq[i+1]-betaseq[i])(log(likePop.Prop)-log(likePop.Acc))))
  
  ## Update rho
  rho <- rho + 1/n*(alpha-0.2338)
  
  ## Update proposed beta
  betaProp <- betaseq[i]/(1+exp(rho))
  
  ## Determine if betaProp has converged
  if(abs(betaProp-betaseq[i+1]) < 0.001){
    betaseq[i+1] <- betaProp
    
    ## Check if all optimal temperatures have been found
    if(betaseq[i+1] > -tBeta*betaMin){
      betaseq[i+1] <- -tBeta*betaMin
      converge <- TRUE
    }
    #If not, propose a new adjacent beta for current converged beta and repeat algorithm
    else{
      i <- i+1
      ## Reinitialize rho,n
      n <- 1
      rho <- initRho
      ## Propose new beta
      betaseq[i+1] <- betaseq[i]/(1+exp(rho))
    }
  }
  else{
    betaseq[i+1] <- betaProp
    n <- n+1
  }
}
