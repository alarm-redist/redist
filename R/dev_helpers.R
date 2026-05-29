# THIS FILE IS NOT INCLUDED IN BUILDS

# call ld_ia() to create a redist_map and redist_plans object for Iowa
# jarl-ignore unused_function: dev-only
ld_ia <- function() {
    data(iowa)
    ia <<- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01)
    plans <<- redist_smc(ia, 100, silent = TRUE)
}