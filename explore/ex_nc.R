library(dplyr)
library(redist)
library(alarmdata)

nc = alarm_50state_map("NC")

#plans = redist_smc(nc, 200, counties=county, runs=4, num_threads_per_process=2, verbose=TRUE)
plans = redist_mergesplit_parallel(nc, 300, counties=county, chains=4, verbose=TRUE)

plans = plans %>%
    mutate(dem = group_frac(nc, ndv, ndv+nrv),
           black = group_frac(nc, vap_black, vap),
           splits = county_splits(nc, county),
           comp = distr_compactness(nc),
           dev = plan_parity(nc))

summary(plans)


