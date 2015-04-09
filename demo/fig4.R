## Replication of Figure 4 in Fifield et. al (2015)
data(algdat.pfull)

######################################################
## 25 precinct, three districts - no pop constraint ##
######################################################
## Get an initial partition
set.seed(1)
initcds <- algdat.pfull$cdmat[,sample(1:ncol(algdat.pfull$cdmat), 1)]

## Run simulations
alg_253 <- redist.mcmc(adjobj = algdat.pfull$adjlist,
                       popvec = algdat.pfull$precinct.data$pop,
                       initcds = initcds,
                       nsims = 10000)

## Get Republican Dissimilarity Index from simulations
rep_dmi_253 <- redist.segcalc(alg_253,
                              algdat.pfull$precinct.data$repvote,
                              algdat.pfull$precinct.data$pop)

## Plot to compare distributions
hist(algdat.pfull$segregation.index$repdiss,
     freq = FALSE,
     cex.lab = 1.2,
     cex.axis = 1.3,
     col = "grey",
     border = "grey",
     xlab = "Republican Dissimilarity Index",
     main = "")
lines(density(rep_dmi_253, from = 0, to = 1))

#######################################################
## 25 precinct, three districts - 20% pop constraint ##
#######################################################
rm(list = ls())
data(algdat.p20)

set.seed(1)
initcds <- algdat.p20$cdmat[,sample(1:ncol(algdat.p20$cdmat), 1)]

## Vector of beta weights
betaweights <- rep(NA, 10); for(i in 1:10){betaweights[i] <- 4^i}

## Run simulations - hard constraint
alg_253_20 <- redist.mcmc(adjobj = algdat.p20$adjlist,
                          popvec = algdat.p20$precinct.data$pop,
                          initcds = initcds,
                          nsims = 10000,
                          popcons = .2)

## Get Republican Dissimilarity Index from simulations - hard constraint
rep_dmi_253_20 <- redist.segcalc(alg_253_20,
                                 algdat.p20$precinct.data$repvote,
                                 algdat.p20$precinct.data$pop)

## Plot to compare distributions
hist(algdat.p20$segregation.index$repdiss,
     freq = FALSE,
     cex.lab = 1.2,
     cex.axis = 1.3,
     col = "grey",
     border = "grey",
     xlab = "Republican Dissimilarity Index",
     main = "")
lines(density(rep_dmi_253_20, from = 0, to = 1))

## Run simulations - tempering population constraint
alg_253_20_st <- redist.mcmc(adjobj = algdat.p20$adjlist,
                             popvec = algdat.p20$precinct.data$pop,
                             initcds = initcds,
                             nsims = 10000,
                             betapop = -5.4,
                             betaweights = betaweights,
                             temperbetapop = 1)

## Resample using inverse probability weighting
alg_253_20_st <- redist.ipw(alg_253_20_st,
                            resampleconstraint = "pop",
                            targetbeta = -5.4,
                            targetpop = .2,
                            temper = 1)

## Get republican dissimilarity index from simulations - soft constraint
rep_dmi_253_20_st <- redist.segcalc(alg_253_20_st,
                                    algdat.p20$precinct.data$repvote,
                                    algdat.p20$precinct.data$pop)

## Add distribution of index using tempering
lines(density(rep_dmi_253_20_st, from = 0, to = 1), col = "red", lty = 5)

#######################################################
## 25 precinct, three districts - 10% pop constraint ##
#######################################################
rm(list = ls())
data(algdat.p10)

set.seed(1)
initcds <- algdat.p10$cdmat[,sample(1:ncol(algdat.p10$cdmat), 1)]

## Vector of beta weights
betaweights <- rep(NA, 10); for(i in 1:10){betaweights[i] <- 4^i}

alg_253_10 <- redist.mcmc(adjobj = algdat.p10$adjlist,
                          popvec = algdat.p10$precinct.data$pop,
                          initcds = initcds,
                          nsims = 10000,
                          popcons = .1)

## Get Republican Dissimilarity Index from simulations
rep_dmi_253_10 <- redist.segcalc(alg_253_10,
                                 algdat.p10$precinct.data$repvote,
                                 algdat.p10$precinct.data$pop)

## Plot to compare distributions
hist(algdat.p10$segregation.index$repdiss,
     freq = FALSE,
     cex.lab = 1.2,
     cex.axis = 1.3,
     col = "grey",
     border = "grey",
     xlab = "Republican Dissimilarity Index",
     main = "")
lines(density(rep_dmi_253_10, from = 0, to = 1))

## Run simulations - tempering population constraint
alg_253_10_st <- redist.mcmc(adjobj = algdat.p10$adjlist,
                             popvec = algdat.p10$precinct.data$pop,
                             initcds = initcds,
                             nsims = 10000,
                             betapop = -9,
                             betaweights = betaweights,
                             temperbetapop = 1)

## Resample using inverse probability weighting
alg_253_10_st <- redist.ipw(alg_253_10_st,
                            resampleconstraint = "pop",
                            targetbeta = -9,
                            targetpop = .1,
                            temper = 1)

## Get republican dissimilarity index from simulations - soft constraint
rep_dmi_253_10_st <- redist.segcalc(alg_253_10_st,
                                    algdat.p10$precinct.data$repvote,
                                    algdat.p10$precinct.data$pop)

## Add distribution of index using tempering
lines(density(rep_dmi_253_10_st, from = 0, to = 1), col = "red", lty = 5)

