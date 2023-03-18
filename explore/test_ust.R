devtools::load_all()

data(iowa)
ia <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01)

bounds = attr(ia, "pop_bounds")
tt1 = sample_ust(ia$adj, ia$pop, bounds[1], bounds[3], rep(1, nrow(ia)))
cat("one\n")

tt2 = sample_ust(ia$adj, ia$pop, bounds[1], bounds[3], as.integer(factor(ia$region)))
cat("two\n")

x = redist_smc(ia, 100, counties=region, silent=TRUE)
cat("map\n")
