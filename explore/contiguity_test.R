devtools::load_all()

# test contigu ----
data(iowa)
iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05)
constr <- redist_constr(iowa_map)
constr <- add_constr_splits(constr, strength = 1.5, admin = name)
constr <- add_constr_grp_hinge(constr, strength = 100,
                               dem_08, tot_08, tgts_group = c(0.5, 0.6))
# encourage districts to have the same number of counties
constr <- constr |>
    add_constr_custom(
        strength = 1000,
        fn = function(plan, distr) {
            # notice that we only use information on precincts in `distr`
            abs(sum(plan == distr) - 99/4)
        }
    )

constr <- constr |>
    add_constr_contiguity(strength = 1000)

redist_mergesplit(map = iowa_map, nsims = 100, warmup = 0, constraints = constr)

redist_flip(map = iowa_map, nsims = 100, warmup = 0, constraints = constr) |>
    dplyr::glimpse()
