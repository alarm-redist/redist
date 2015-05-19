context("Test functionality of redist.rsg")

rm(list = ls())

## Get set of initial cds
data(algdat.pfull)
set.seed(1)
cds1 <- algdat.pfull$cdmat[,sample(1:ncol(algdat.pfull$cdmat), 1)]

## Set parameters
ndists <- 3
targpop <- sum(algdat.pfull$precinct.data$pop) / ndists

test_that("redist.rsg works", {
    
    ## Test basic redist.rsg
    test1 <- redist.rsg(adj.list = algdat.pfull$adjlist,
                        population = algdat.pfull$precinct.data$pop,
                        ndists = ndists,
                        thresh = .2)
    expect_equal(length(unique(test1$district_membership)), ndists)
    
    ## Test in combination with redist.mcmc
    test2 <- redist.mcmc(adjobj = algdat.pfull$adjlist,
                         popvec = algdat.pfull$precinct.data$pop,
                         nsims = 100,
                         ndists = ndists)
    expect_equal(length(unique(test2$partitions[100])), ndists)

})

