devtools::load_all(here::here("."))
library(tictoc)

data(iowa)
ia = redist_map(iowa, existing_plan=cd_2010, pop_tol=0.01)

cat("SINGLE-THREADED:\n")
tic()
plans = redist_smc(ia, 10000, cores=1, verbose=F)
toc()

cat("\n\nMULTIPLE THREADS:\n")
tic()
plans = redist_smc(ia, 10000, cores=4, verbose=F)
toc()

summary(plans)
