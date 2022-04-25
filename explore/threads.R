devtools::load_all(here::here("."))
library(tictoc)

data(iowa)
ia = redist_map(iowa, existing_plan=cd_2010, pop_tol=0.01)

cat("SINGLE-THREADED:\n")
tic()
plans = redist_smc(ia, 5000, cores=1, verbose=T)
toc()

cat("\n\nMULTIPLE THREADS:\n")
tic()
plans = redist_smc(ia, 5000, cores=2, verbose=T)
toc()

summary(plans)
