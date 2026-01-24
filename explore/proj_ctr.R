library(alarmdata)
library(sf)
library(tidyverse)
library(ggredist)
devtools::load_all("../ggredist/")
devtools::load_all(".")

centroid <- function(plans, map) {
    X <- st_centroid(map) |>
        suppressWarnings() |>
        st_coordinates() |>
        cbind(pop=map$pop) |>
        as.data.frame()

    out <- matrix(nrow(plans), ncol=2)
    m <- as.matrix(plans)
    nd <- attr(plans, "ndists")
    fun <- \(x) c(weighted.mean(x$X, x$pop), weighted.mean(x$Y, x$pop))
    apply(m, 2, function(pl) {
        tapply(X, pl, fun) |>
            do.call(rbind, args=_)
    }, simplify=FALSE) |>
        do.call(rbind, args=_) |>
        `colnames<-`(c("ctr_x", "ctr_y")) |>
        as_tibble()
}
# Fisher's method to combine p-values
pv_combine <- function(...) {
    pv = do.call(cbind, list(...))
    pchisq(-2 * rowSums(log(pv)), df = 2*ncol(pv))
}

d_ccd = tigris::county_subdivisions("SC", year=2020, cb=TRUE)
map = alarm_50state_map("SC", year=2020) |>
    mutate(ccd = d_ccd$NAME[geomander::geo_match(map, d_ccd, "center")])
plans = alarm_50state_plans("SC", year=2020) |>
    mutate(centroid(pl(), map))

plot(map, proj_avg(plans, vap_black/total_vap)) +
    geom_district(aes(group = cd_2020), fill=NA, is_coverage=TRUE) +
    crop_to(county %in% c("101", "045", "029", "091", "017"), map) +
    # scale_fill_party_c(limits=c(0.3, 0.7)) +
    theme(legend.key.width=unit(1, "cm"))

pc_x = proj_contr(plans, ctr_x, pfdr=FALSE)
pc_y = proj_contr(plans, ctr_y, pfdr=FALSE)
pc_pv = proj_contr(plans, ctr_x^2 + ctr_y^2, pfdr=TRUE)
pc_dem = proj_contr(plans, ndshare, pfdr=TRUE)
pc_bvap = proj_contr(plans, vap_black/total_vap, pfdr=TRUE)
ctr = st_centroid(map) |>
    suppressWarnings() |>
    st_coordinates() |>
    as_tibble() |>
    mutate(
        pc_x = pc_x,
        pc_y = pc_y,
        pa_x = proj_avg(plans, ctr_x),
        pa_y = proj_avg(plans, ctr_y),
        q_disp = attr(pc_pv, "q"),
        # q_dem = attr(pc_dem, "q"),
        q_bvap = attr(pc_bvap, "q"),
    )

redist.plot.contr_pfdr(map, pc_bvap, level=0.1) +
    geom_district(aes(group = cd_2020), fill=NA, is_coverage=TRUE) +
    # crop_to(county %in% c("101", "045", "029", "091", "017"), map) +
    # scale_fill_party_c(name="Vote share vs. simulations\n",
    #                              midpoint=0, limits=c(-0.2, 0.2),
    #                              labels=ggredist::label_party_margin(midpoint=0)) +
    wacolors::scale_fill_wa_c("vantage", name="BVAP vs. simulations",
                              midpoint=0, limits=c(-0.25, 0.25)) +
    theme(legend.key.width=unit(1.5, "cm"))

redist.plot.contr_pfdr(map, pc_pv, level=0.5) +
    wacolors::scale_fill_wa_c("vantage", midpoint=0)

# map$geometry = sf::st_make_valid(map$geometry)
map |>
    bind_cols(ctr) |>
    group_by(cd_2020, ccd) |>
    st_as_sf() |>
    summarize(across(X:pa_y, weighted.mean, w=vap),
              across(q_disp:q_bvap, ~ pchisq(-2*sum(log(.)), 2*n(), lower.tail=FALSE)),
              across(c(starts_with("vap"), area_land), sum),
              is_coverage=TRUE, .groups="drop") |>
ggplot() +
    # geom_sf(aes(fill = ndv / (ndv + nrv)), linewidth=0, color=NA) +
    geom_sf(aes(fill = vap_black/vap, alpha=1609^2*vap/area_land), linewidth=0, color=NA) +
    geom_district(aes(group = cd_2020), fill=NA, is_coverage=TRUE) +
    # geom_segment(aes(x=X, y=Y, xend=X+pc_x/4, yend=Y+pc_y/4),
    # geom_segment(aes(xend=X, yend=Y, x=X-pc_x/5, y=Y-pc_y/5),
    geom_segment(aes(x=X, y=Y, xend=X+pc_x/5, yend=Y+pc_y/5),
    # geom_segment(aes(x=X, y=Y, xend=X + (pa_x-X)/5, yend=Y + (pa_y-Y)/5),
                 # data = filter(ctr, q_disp < 0.4),
                 # data = filter(ctr, q_bvap < 0.05),
                 data = \(x) filter(x, q_bvap < 0.01),
                 # data = ctr,
                 # linewidth=0.25, alpha=0.5,
                 linewidth=0.4, alpha=0.5,
                 arrow=arrow(angle=30, length=unit(3, "pt"), type="closed")) +
    # crop_to(county %in% c("101", "045", "029", "091", "017"), map) +
    # crop_to(county %in% c("101"), map) +
    # scale_fill_party_c() +
    # stat_cities(geom="text", color="white", adjust=0.8) +
    scale_alpha_continuous(range=c(0.1, 1.0), trans="sqrt", limits=c(0, 2000),
                           oob=scales::squish, guide="none") +
    wacolors::scale_fill_wa_c("vantage", name="BVAP",
                              midpoint=0.5, limits=c(0.2, 0.8), oob=scales::squish) +
    # wacolors::scale_fill_wa_c("volcano", reverse=T) +
    theme_void() +
    theme(legend.position="bottom",
          legend.key.width=unit(2, "cm"))


# idx_nola = which(map$county == "Orleans Parish")
# overl <- function(plans, idx) {
#     nd = attr(plans, "ndists")
#     apply(as.matrix(plans), 2, function(pl) {
#         tabulate(pl[idx], nbins=nd)/tabulate(pl, nbins=nd)
#     }) |>
#         as.numeric()
# }
# plans$overl_nola = overl(plans, idx_nola)
#
# plot(plans, overl_nola, geom="boxplot")
# redist.plot.contr_pfdr(map, proj_contr(plans, overl_nola, pfdr=TRUE), level=0.01) +
#     wacolors::scale_fill_wa_c("vantage", midpoint=0)
#
# subset_sampled(plans) |>
#     plan_distances(ncores=4) |>
#     classify_plans() |>
#     plot()

co = prec_cooccurrence(plans, ncores=4)
m_cons = lapply(seq_len(attr(plans, "ndists")), function(i) {
    # matrixStats::colWeightedMeans(co[map$cd_2010 == i, ], map$vap[map$cd_2010 == i])
    matrixStats::rowMeans2(as.matrix(plans) == i)
}) |>
    do.call(cbind, args=_)
pr_cons = m_cons[cbind(seq_len(nrow(map)), max.col(m_cons))]

ggplot(map, aes(fill=as.factor(max.col(m_cons)), alpha=pr_cons)) +
    geom_sf(linewidth=0, color=NA) +
    geom_district(aes(group = cd_2020), alpha=1, fill=NA, is_coverage=TRUE) +
    scale_alpha_continuous(range=c(0.0, 1.0), limits=c(0, 1),
                           oob=scales::squish, guide="none") +
    wacolors::scale_fill_wa_d("skagit", guide="none") +
    theme_void()


edges = dplyr::as_tibble(map) %>%
    sf::st_as_sf() %>%
    dplyr::select(geometry = attr(map, "sf_column")) %>%
    sf::st_intersection() %>%
    dplyr::as_tibble() %>%
    dplyr::filter(.data$n.overlaps == 2) %>%
    dplyr::mutate(i = sapply(.data$origins, function(x) x[1]),
                  j = sapply(.data$origins, function(x) x[2])) %>%
    dplyr::filter(sf::st_dimension(.data$geometry) == 1) %>%
    dplyr::select(i, j, geometry) |>
    sf::st_as_sf()

edge_cuts = imap(map$adj, function(nbors, i) {
    new_tibble(list(
        i = rep(i, length(nbors)),
        j = nbors + 1,
        cut = 1 - co[i, nbors + 1]
    ))
}) |>
    bind_rows() |>
    filter(i < j) |>
    left_join(edges, by=c("i", "j")) |>
    st_as_sf()

ggplot(map) +
    # geom_district(aes(group=1), is_coverage=TRUE, fill=NA) +
    # geom_sf(aes(fill=vap_black/vap, color=vap_black/vap), linewidth=0.1) +
    geom_sf(aes(linewidth=cut, alpha=cut), data=filter(edge_cuts, cut > 0.05),
            color="black") +
    geom_district(aes(group=cd_2010), color="red", linewidth=0.1,
                  is_coverage=TRUE, fill=NA) +
    scale_linewidth_continuous(range=c(0.0, 1.1), limits=0:1, guide="none") +
    scale_alpha_continuous(range=c(0.8, 1), limits=0:1, guide="none") +
    wacolors::scale_fill_wa_c("vantage", name="BVAP", labels=scales::percent,
                              midpoint=0.5, limits=0:1, which=2:14) +
    wacolors::scale_color_wa_c("vantage", name="BVAP", labels=scales::percent,
                              midpoint=0.5, limits=0:1, which=2:14) +
    # annotate("text", x=I(0.15), y=I(0.25),
    #          label="White = Pr(edge cut)\nBlack=2020 districts") +
    theme_void() +
    theme(legend.position=c(0.85, 0.15))


X <- st_centroid(map) |>
    suppressWarnings() |>
    st_coordinates()
dm = st_distance(st_centroid(map))
m <- as.matrix(plans)[, -1]
N = ncol(m)
disp = matrix(0, nrow=nrow(map), ncol=2)
for (pl in 1:ncol(m)) {
    for (i in 1:nrow(m)) {
        j = which.min(dm[i, m[, pl] != m[i, pl]])
        disp[i, ] = disp[i, ] + (X[j, ] - X[i, ])/N
    }
}
colnames(disp) = c("dX", "dY")

ggplot(map) +
    geom_sf(aes(fill = vap_black/vap, alpha=1609^2*vap/area_land), linewidth=0, color=NA) +
    geom_district(aes(group = cd_2020), fill=NA, is_coverage=TRUE) +
    geom_segment(aes(x=X, y=Y, xend=X-dX/7, yend=Y-dY/7),
                 data = as.data.frame(cbind(X, disp)),
                 linewidth=0.4, alpha=0.5,
                 arrow=arrow(angle=30, length=unit(3, "pt"), type="closed")) +
    theme_void()
