# Understand how VI scales -----
library(sf)


make_map = function(size=4, n_distr=2) {
    geom = st_sfc(st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,0))))) %>%
        st_make_grid(n=c(size, size)) %>%
        st_set_crs(4269)
    d = st_as_sf(geom, pop=1)
    redist_map(d, ndists=n_distr, pop_tol=0.1)
}

params = crossing(size=c(4, 8, 16, 32), n_distr=c(2, 4, 8, 16, 32, 64)) %>%
    filter(size^2/n_distr >= 2)

res = pmap_dfr(params, function(size, n_distr) {
    map = make_map(size, n_distr)
    plans = redist_smc(map, 500, silent=TRUE)
    div = plans_diversity(plans)
    tibble(size = size, n_distr = n_distr,
           min = min(div),
           max = max(div),
           mean = mean(div))
})

ggplot(res, aes(n_distr, max, color=size, group=size)) +
    geom_function(fun=log, inherit.aes=F, color="gray") +
    #geom_function(fun=\(x) log(x)-log(size), inherit.aes=F, color="gray") +
    geom_line(aes(x=n_distr, y=log(n_distr)), color="gray") +
    geom_line() +
    geom_point() +
    scale_x_continuous(trans="log2") +
    scale_y_continuous(trans="log2") +
    scale_color_continuous(trans="log2")


map = make_map(16, 16)
plans = redist_smc(map, 5000)
div = plans_diversity(plans, n_max=500)
tibble(size = size, n_distr = n_distr,
       min = min(div),
       max = max(div),
       mean = mean(div))

