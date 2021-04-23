devtools::load_all(".")

data(fl25)
data(iowa)
adj = list(1L:14L, c(0L, 3L, 15L, 16L, 17L, 18L, 19L, 20L, 21L), c(0L, 3L, 5L, 20L, 22L),
           c(0L, 1L, 2L, 20L), c(0L, 6L, 21L, 23L, 24L), c(0L, 2L, 8L, 22L), c(0L, 4L, 7L),
           c(0L, 6L, 9L), c(0L, 5L, 13L), c(0L, 7L, 10L), c(0L, 9L, 11L), c(0L, 10L, 12L),
           c(0L, 11L, 14L), c(0L, 8L, 14L), c(0L, 12L, 13L), c(1L, 16L, 17L), c(1L, 15L, 18L),
           c(1L, 15L, 20L), c(1L, 16L, 19L), c(1L, 18L, 21L), c(1L, 2L, 3L, 17L, 22L),
           c(1L, 4L, 19L, 23L, 24L), c(2L, 5L, 20L), c(4L, 21L, 24L), c(4L, 21L, 23))
pop = fl25$pop

if (F) {
#N = 20
#res = redist.smc(adj, pop, N, 3, pop_tol=0.1, silent=T)

set.seed(5118)
ia_adj = redist.adjacency(iowa)
ust = sample_ust(ia_adj, iowa$pop, 723509, 799668, as.integer(as.factor(iowa$region)))

{
centers = suppressWarnings(st_centroid(iowa))
nb <- lapply(ust, function(x) x + 1L)
edgedf <- tibble(start = rep(1:length(nb), lengths(nb)),
                 finish = unlist(nb)) %>%
    rowwise() %>%
    mutate(geometry = st_sfc(st_linestring(matrix(
        c(
            as.numeric(centers$geometry[[start]]),
            as.numeric(centers$geometry[[finish]])
        ),
        nrow = 2,
        byrow = TRUE
    ))))
edge_coord = tibble::as_tibble(st_coordinates(edgedf$geometry)) %>%
    mutate(which = rep(c("start", "end"), n()/2)) %>%
    tidyr::pivot_wider(names_from=which, values_from=c(X, Y))
edgedf = bind_cols(edgedf, edge_coord)


ggplot(iowa, aes(fill = region)) +
    geom_sf() +
    geom_sf(data=centers, inherit.aes=F) +
    ggplot2::geom_segment(aes(x=X_start, y=Y_start, xend=X_end, yend=Y_end),
                          arrow = grid::arrow(angle=25, length=grid::unit(0.1, 'inches'), type='closed'),
                          data=edgedf, inherit.aes=F) +
    theme_void()
}
}

iowa_map = redist_map(iowa, ndists=4, pop_tol=0.05)
plans = redist_smc(iowa_map, 10000, counties=region) %>%
    mutate(dem = group_frac(iowa_map, dem_08, tot_08))
plans = redist_smc(iowa_map, 10000) %>%
    mutate(dem = group_frac(iowa_map, dem_08, tot_08))

pa_shp = sf::read_sf("../smc-paper/data/PA/PA/pa_final.shp")
pa = redist_map(pa_shp, ndists=18, total_pop=POP100, pop_tol=0.005)

{
tic()
pl = redist_smc(pa, 100, counties=COUNTYFP10)
toc()
}
