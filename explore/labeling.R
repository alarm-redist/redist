library(ggplot2)

shp = geomander::checkerboard %>%
    filter(i<4, j<4) %>%
    mutate(pop=1) %>%
    sf::st_set_crs(3857)
map = redist_map(shp, ndists=4, pop_tol=0.1)

to_lbl = function(x) paste(x, collapse="")

plans = redist_smc(map, 100e3, compactness=0, seq_alpha=0.5, adjust_labels=T, truncate=F, trunc_fn=NULL)

renumb = apply(as.matrix(plans), 2, \(x) order(unique(x)))
m_renumb = redist:::renumber_matrix(as.matrix(plans), renumb)

smc_unlab = as.integer(as.factor(apply(m_renumb, 2, to_lbl)))
plot(table(smc_unlab))

nlab = tibble(pl=smc_unlab, nlab=exp(attr(plans, "diagnostics")$log_labels)) %>%
    group_by(pl) %>%
    summarize(nlab = mean(nlab)) %>%
    dplyr::pull(nlab)

log_st = tibble(pl=smc_unlab, log_st=by_plan(comp_log_st(plans, map))) %>%
    group_by(pl) %>%
    summarize(log_st = mean(log_st)) %>%
    dplyr::pull(log_st)

qplot(table(smc_unlab), 2*nlab)
qplot(x1, table(smc_unlab), color=nlab)
qplot(x1, table(smc_unlab), color=sqrt(log_st))
qplot(log(x3)-log_st, log(x4)-log_st, color=sqrt(log_st))
qplot(log(x3)-log_st, log(x6)-log_st, color=sqrt(log_st))
#log_st[which(!duplicated(smc_unlab))]

which.max(x1)
redist.plot.plans(plans, which.max(smc_unlab==64), shp)
