devtools::load_all(".", recompile=F)

data(algdat.p10)

m = redist.smc(algdat.p10$adjlist, algdat.p10$precinct.data$pop,
               nsims=100, ndists=3, popcons=0.1, compactness=0.5)
print(m)
cat('\n')

m = redist.smc(algdat.p10$adjlist, algdat.p10$precinct.data$pop,
               nsims=100, ndists=3, popcons=0.1, compactness=0.0)
print(m)
