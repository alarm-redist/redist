context("Test functionality of redist.mcmc")

rm(list = ls())

## Get set of initial cds
data(algdat.pfull)
set.seed(1)
cds1 <- algdat.pfull$cdmat[,sample(1:ncol(algdat.pfull$cdmat), 1)]
nsims <- 100

## Initial tests
test_that("redist.mcmc works", {

    ## Basic simulator
    test1 <- redist.mcmc(adjobj = algdat.pfull$adjlist,
                         popvec = algdat.pfull$precinct.data$pop,
                         initcds = cds1,
                         nsims = nsims)
    
    expect_is(test1, "redist")

    expect_equal(length(test1), 9)
    expect_equal(ncol(test1$partitions), nsims)
    expect_equal(length(test1$distance_parity), nsims)
    expect_equal(length(test1$mhdecisions), nsims)
    expect_equal(length(test1$mhprob), nsims)
    expect_equal(length(test1$pparam), nsims)
    expect_equal(length(test1$constraint_pop), nsims)
    expect_equal(length(test1$constraint_compact), nsims)
    expect_equal(length(test1$constraint_segregation), nsims)
    expect_equal(length(test1$constraint_similar), nsims)

    ## Population constraint
    test2 <- redist.mcmc(adjobj = algdat.pfull$adjlist,
                         popvec = algdat.pfull$precinct.data$pop,
                         initcds = cds1,
                         nsims = nsims,
                         popcons = .2)

    expect_less_than(max(test2$distance_parity), max(test1$distance_parity))

    ## With soft population constraint
    test3 <- redist.mcmc(adjobj = algdat.pfull$adjlist,
                         popvec = algdat.pfull$precinct.data$pop,
                         initcds = cds1,
                         nsims = nsims,
                         betapop = -5)

    expect_less_than(max(test3$distance_parity), max(test1$distance_parity))

    ## Soft population constraint and tempering
    test4 <- redist.mcmc(adjobj = algdat.pfull$adjlist,
                         popvec = algdat.pfull$precinct.data$pop,
                         initcds = cds1,
                         nsims = nsims,
                         betapop = -5,
                         temperbetapop = 1)

    expect_equal(length(test4), 12)
    expect_equal(length(test4$beta_sequence), nsims)
    expect_equal(length(test4$mhdecisions_beta), nsims)
    expect_equal(length(test4$mhprob_beta), nsims)

    expect_more_than(length(unique(test4$beta_sequence)), 1)

    ## Soft compactness constraint and tempering
    test5 <- redist.mcmc(adjobj = algdat.pfull$adjlist,
                         popvec = algdat.pfull$precinct.data$pop,
                         initcds = cds1,
                         nsims = nsims,
                         ssdmat = algdat.pfull$distancemat,
                         betacompact = -2,
                         temperbetacompact = 1)

    expect_more_than(length(unique(test5$beta_sequence)), 1)

    ## Soft segregation constraint and tempering
    test6 <- redist.mcmc(adjobj = algdat.pfull$adjlist,
                         popvec = algdat.pfull$precinct.data$pop,
                         initcds = cds1,
                         grouppopvec = algdat.pfull$precinct.data$repvote,
                         nsims = nsims,
                         betaseg = -2,
                         temperbetaseg = 1)

    expect_more_than(length(unique(test6$beta_sequence)), 1)

    ## Soft similarity constraint and tempering
    test7 <- redist.mcmc(adjobj = algdat.pfull$adjlist,
                         popvec = algdat.pfull$precinct.data$pop,
                         initcds = cds1,
                         grouppopvec = algdat.pfull$precinct.data$repvote,
                         nsims = nsims,
                         betasimilar = -2,
                         temperbetasimilar = 1)

    expect_more_than(length(unique(test7$beta_sequence)), 1)

    ## Multiple swaps test
    test8 <- redist.mcmc(adjobj = algdat.pfull$adjlist,
                         popvec = algdat.pfull$precinct.data$pop,
                         initcds = cds1,
                         lambda = 1,
                         nsims = nsims)
    expect_more_than(max(test8$pparam), 1)

    ## Multiple save points test
    test9 <- redist.mcmc(adjobj = algdat.pfull$adjlist,
                         popvec = algdat.pfull$precinct.data$pop,
                         initcds = cds1,
                         lambda = 1,
                         nsims = nsims,
                         nloop = 2,
                         savename = "test")
    load("test.RData")

    expect_is(algout, "redist")
    expect_equal(ncol(algout$partitions), nsims * 2)

    ## Algorithm restart test
    test10 <- redist.mcmc(adjobj = algdat.pfull$adjlist,
                          popvec = algdat.pfull$precinct.data$pop,
                          initcds = cds1,
                          lambda = 1,
                          nsims = nsims,
                          nloop = 4,
                          loopscompleted = 2,
                          savename = "test")
    load("test.RData")

    expect_is(algout, "redist")
    expect_equal(ncol(algout$partitions), nsims * 4)

    ## Test rng functionality
    set.seed(1)
    test11a <- redist.mcmc(adjobj = algdat.pfull$adjlist,
                           popvec = algdat.pfull$precinct.data$pop,
                           initcds = cds1,
                           nsims = nsims)
    
    set.seed(1)
    test11b <- redist.mcmc(adjobj = algdat.pfull$adjlist,
                           popvec = algdat.pfull$precinct.data$pop,
                           initcds = cds1,
                           nsims = nsims)

    expect_identical(test11a, test11b)

    ## Test rng functionality with multiple loops
    set.seed(1)
    test12a <- redist.mcmc(adjobj = algdat.pfull$adjlist,
                           popvec = algdat.pfull$precinct.data$pop,
                           initcds = cds1,
                           nloop = 2,
                           nsims = nsims,
                           savename = "testmsa")
    
    set.seed(1)
    test12b <- redist.mcmc(adjobj = algdat.pfull$adjlist,
                           popvec = algdat.pfull$precinct.data$pop,
                           initcds = cds1,
                           nloop = 2,
                           nsims = nsims,
                           savename = "testmsb")
    load("testmsa.RData")
    msa <- algout
    load("testmsb.RData")
    msb <- algout

    expect_identical(msa, msb)
    
    
})
