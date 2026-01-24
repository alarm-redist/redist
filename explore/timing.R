devtools::install(quick=T, upgrade=F)
devtools::load_all(".")

data(iowa)
ia = redist_map(iowa, existing_plan=cd_2010, pop_tol=0.001)

bench::mark(
    SMC = redist_smc(ia, 200, silent=TRUE),
    MS = redist_mergesplit(ia, 200, silent=TRUE),
    check = FALSE,
    min_iterations = 10,
    memory = FALSE
)
