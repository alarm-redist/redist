devtools::load_all(here::here("."))
library(tictoc)

data(iowa)
ia = redist_map(iowa, existing_plan=cd_2010, pop_tol=0.01)

cat("SINGLE-THREADED:\n")
tic()
plans = redist_smc(ia, 5000, ncores=1, verbose=T)
toc()

cat("\n\nMULTIPLE THREADS:\n")
tic()
plans = redist_smc(ia, 5000, ncores=2, verbose=F)
toc()

summary(plans)
