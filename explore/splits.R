devtools::load_all(".")
library(patchwork)

data(nh_map, package="redistmetrics")

nh_map$muni <- substr(nh_map$vtd, 1, 4)

cons <- redist_constr(nh_map) %>%
    add_constr_splits(strength = 20, admin = muni)

set.seed(12345)
pl_no_con <- redist_mergesplit(nh_map, nsims = 200, init_plan = nh_map$r_2020) %>%
    mutate(spl = splits_admin(pl(), nh_map, muni))
set.seed(12345)
pl_con <- redist_mergesplit(nh_map, nsims = 200, init_plan = nh_map$r_2020,
                            constraints = cons) %>%
    mutate(spl = splits_admin(pl(), nh_map, muni))

(pl_no_con %>%
    redist.plot.trace(spl) +
    labs(title = 'MS without constraint')) +
(pl_con %>%
    redist.plot.trace(spl) +
    labs(title = 'MS with constraint') +
    plot_annotation(title = 'ba7c109'))



set.seed(12345)
pl_no_con <- redist_smc(nh_map, nsims = 200) %>%
    mutate(spl = splits_admin(pl(), nh_map, muni))
set.seed(12345)
pl_con <- redist_smc(nh_map, nsims = 200, constraints = cons) %>%
    mutate(spl = splits_admin(pl(), nh_map, muni))

(pl_no_con %>%
        redist.plot.trace(spl) +
        labs(title = 'SMC without constraint')) +
(pl_con %>%
     redist.plot.trace(spl) +
     labs(title = 'SMC with constraint') +
     plot_annotation(title = 'ba7c109'))
