#################################
# Author: Ben Fifield
# Created: 2014/09/21
# Last Revision: 2014/09/22
# Institution: Princeton University
# Purpose: Called by runSWA.R after combineDataSplit.R
#          to generate diagnostics
#################################

## Number of acceptances
accept <- 0
for(z in 2:length(edgecuts[[10]])){
    if(edgecuts[[10]][z-1] != edgecuts[[10]][z]){
        accept <- accept + 1
    }
}

## Distance from population parity
paritydist <- distParity(edgecuts[[1]], pcData$pop, parity)

## For subsetting accept/reject plot
ind <- seq(1, length(edgecuts[[16]]), by = 25)

######################
## Diagnostic plots ##
######################
pdf(file = paste(dwd, "diagPlot", state, "_", (1 - eprob) * 100,
        "_", margin.pct, "_", lambda,
        "_b", params$initbeta[aid] * -1, "a", abeta,
        "_bDiss", params$initbetadiss[aid], "a", abetadiss,
        "_bPop", params$initbetapop[aid] * -1, "a", abetapop,
        "t", params$targbetapop[aid] * -1, "b", params$bybetapop[aid] * -1,
        "_bSwitch", params$initbetaswitch[aid] * -1, "a", abetaswitch,
        "_pow", wpow,
        "_par", pnum, ".pdf", sep = ""),
    height = 6, width = 12)
par(mfrow = c(1,2))
## Trace of beta
plot(edgecuts[[10]], type = "l",
     xlab = "Iterations", ylab = "",
     main = paste("Trace of Beta \n Lambda =",
         params$lambda[aid], "pow =",
         params$weightpow[aid], sep = " "))
## Distribution of beta values
barplot(table(edgecuts[[10]]),
        main = paste("Distribution of Beta_pop \n Acceptance Rate =",
            accept / length(edgecuts[[10]]), sep = " "))

## Accept/Reject
plot(edgecuts[[16]][ind], pch = 16,
     cex = .5,
     yaxt = "n",
     col = ifelse(edgecuts[[16]][ind] == 0, "red", "black"),
     xlab = "Iterations", ylab = "",
     main = paste("Rejections from SWA \n Acceptance Rate =",
         sum(edgecuts[[16]]) / length(edgecuts[[16]]), sep = " "))
axis(2, c(0, 1), c("Reject", "Accept"))
## Trace of district population
plot(paritydist, type = "l",
     xlab = "Iterations", ylab = "",
     main = "Trace of Distance from Parity")
dev.off()



