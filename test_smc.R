devtools::load_all(".", recompile=F)

data(algdat.p10)

gctorture(TRUE)
cat("TEST 1\n")
redist.smc(algdat.p10$adjlist, algdat.p10$precinct.data$pop,
           nsims=100, ndists=3, popcons=0.1, compactness=1)
gc()


cat("TEST 2B\n")
redist.smc(algdat.p10$adjlist, algdat.p10$precinct.data$pop,
           nsims=1000, ndists=3, popcons=0.1, compactness=1)
gc()

cat("TEST 3\n")
redist.smc(algdat.p10$adjlist, algdat.p10$precinct.data$pop,
           nsims=10000, ndists=3, popcons=0.1, compactness=1)
gc()

cat("TEST 4\n")
redist.smc(algdat.p10$adjlist, algdat.p10$precinct.data$pop,
           nsims=100, ndists=3, popcons=0.1, compactness=0.5)
gc()

redist.smc(algdat.p10$adjlist, algdat.p10$precinct.data$pop,
           nsims=100, ndists=3, popcons=0.1, compactness=0.0)
gc()
