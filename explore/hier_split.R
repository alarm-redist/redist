devtools::load_all(".")
library(purrr)

ld_ia()
V = nrow(ia)

all_below <- function(i, ust) {
    children = ust[[i]]
    if (length(children) == 0)
        i
    else
        c(i, unlist(map(children + 1, all_below, ust)))
}

gen_plan <- function(...) {
    plan = rep(0, V)

    # split 1
    tgt_1 = get_target(ia) * 2
    ust = sample_ust(ia$adj, ia$pop, get_target(ia), get_target(ia), rep(1, V), rep(F, V))
    pop_below = map_dbl(seq_along(ust), ~ tree_pop(ust, . - 1, ia$pop, rep(0, V), rep(0, V)))

    idx = all_below(which.min(abs(pop_below - tgt_1)), ust)
    plan[idx] = 1
    plan[-idx] = 3
    # plot(ia, plan)

    # split 2
    tgt_2 = get_target(ia)
    sub = plan == 1
    ust = sample_ust(ia$adj, ia$pop, get_target(ia), get_target(ia), rep(1, V), !sub)
    sub = which(sub)
    pop_below = map_dbl(sub, ~ tree_pop(ust, . - 1, ia$pop, rep(0, V), rep(0, V)))
    idx = all_below(sub[which.min(abs(pop_below - tgt_2))], ust)
    plan[idx] = 1
    plan[setdiff(sub, idx)] = 2
    # plot(ia, plan)

    # split 3
    sub = plan == 3
    ust = sample_ust(ia$adj, ia$pop, get_target(ia), get_target(ia), rep(1, V), !sub)
    sub = which(sub)
    pop_below = map_dbl(sub, ~ tree_pop(ust, . - 1, ia$pop, rep(0, V), rep(0, V)))
    idx = all_below(sub[which.min(abs(pop_below - tgt_2))], ust)
    plan[idx] = 3
    plan[setdiff(sub, idx)] = 4
    # plot(ia, plan)

    max(abs(tapply(ia$pop, plan, sum) / get_target(ia) - 1))
}

devs = map_dbl(1:1000, gen_plan, .progress=TRUE)

map = alarmdata::alarm_50state_map("wa")
V = nrow(map)
bounds = attr(map, "pop_bounds") %o% 1:9
rev_bounds = sum(map$pop) - bounds
bounds = rbind(
    pmax(bounds[1, ], rev(rev_bounds[3, ])),
    bounds[2, ],
    pmin(bounds[3, ], rev(rev_bounds[1, ]))
)
mean(map$pop < 0.01*get_target(map))

tgt_1 = sum(map$pop)/2
ust = sample_ust(map$adj, map$pop, tgt_1, tgt_1, vctrs::vec_group_id(map$county), rep(F, V))
ust = sample_ust(map$adj, map$pop, tgt_1*0.1, tgt_1*1.9, vctrs::vec_group_id(map$county), rep(F, V))
# ust = sample_ust(map$adj, map$pop, tgt_1, tgt_1, rep(1, V), rep(F, V))
pop_below = map_dbl(seq_along(ust), ~ tree_pop(ust, . - 1, map$pop, rep(0, V), rep(0, V)))
plot(sort(pop_below), cex=0.3); abline(h=bounds, col='red')
