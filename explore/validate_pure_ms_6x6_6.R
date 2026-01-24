# NOTE: N was changed to M
library(here)
devtools::load_all(recompile=F)
# library(redist)
library(sf)
library(tidyverse)
library(patchwork)
library(scales)
# source(here("validation/R/setup.R"))

PAL = c("#CC79A7", "#E69F00", "#56B4E9", "#20B073", "#0072B2", "#D55E00", "#999999")

n = 6

shp = geomander::checkerboard |>
    filter(i<n, j<n) |>
    st_set_crs(3857) |>
    mutate(pop = n)
map = redist_map(shp, ndists=n, pop_tol=0.01)


# enumerate all valid plans -----
path_enum <- here(sprintf("../smc-paper/julia/trees/enumerations/%dx%d_%d.txt", n,n,n))
enum_raw = read_csv(path_enum, col_names=FALSE, col_types=cols(.default="i")) %>%
    as.matrix() %>%
    t()
enum_raw = list(plans=enum_raw)
log_st_enum = by_plan(comp_log_st(enum_raw$plans, map))
w = exp(log_st_enum - max(log_st_enum))
w = w / sum(w)

pl_enum = redist_plans(enum_raw$plans, map, "enumpart") %>%
    mutate(edges = comp_edges_rem(pl(), map))

true_median = weighted.mean(by_plan(pl_enum$edges == 22, n), w)
true_htail = weighted.mean(by_plan(pl_enum$edges >= 27, n), w)
true_ltail = weighted.mean(by_plan(pl_enum$edges == 18, n), w)

# Run comparisons --------
M_seq = 10 * round(10^seq(0, 3, by=0.5))

pls_smc = list()

set.seed(5118)
res_smc = imap_dfr(M_seq, function(M, i) {
    cat("Running SMC for M =", M, "\n")

    pl_smc = map(1:20, function(j) {
        redist_mergesplit(map, 100 + M*10, warmup=100, thin=10, silent=TRUE) |>
            mutate(edges = comp_edges_rem(pl(), map))
    }) |>
        do.call(rbind, args=_)
    pls_smc[[i]] <<- pl_smc # save it

    diagn = attr(pl_smc, "diagnostics")
    vi = plans_diversity(pl_smc, 200L)
    chain = by_plan(pl_smc$chain, ndists=n)
    edges = by_plan(pl_smc$edges, ndists=n)

    ests_median = tapply(edges == 22, chain, mean)
    ests_htail = tapply(edges >= 27, chain, mean)
    ests_ltail = tapply(edges == 18, chain, mean)

    tibble_row(alg="smc", iter=M,
               # time = mean(map_dbl(diagn, ~ .$runtime)),
               time = mean(unlist(diagn)),
               # n_eff = mean(map_dbl(diagn, ~ .$n_eff)),
               vi_q10 = quantile(vi, 0.1),
               vi_q90 = quantile(vi, 0.9),
               rhat_median = redist:::diag_rhat(edges == 22, chain),
               rhat_htail = redist:::diag_rhat(edges >= 27, chain),
               rhat_ltail = redist:::diag_rhat(edges == 18, chain),
               sd_median = sd(ests_median),
               sd_htail = sd(ests_htail),
               sd_ltail = sd(ests_ltail),
               est_median = mean(ests_median),
               est_htail = mean(ests_htail),
               est_ltail = mean(ests_ltail),
               bias_median = est_median - true_median,
               bias_htail = est_htail - true_htail,
               bias_ltail = est_ltail - true_ltail,
               rmse_median = sqrt(mean((ests_median - true_median)^2)),
               rmse_htail = sqrt(mean((ests_htail - true_htail)^2)),
               rmse_ltail = sqrt(mean((ests_ltail - true_ltail)^2)))
})
names(pls_smc) = as.character(M_seq)

res = res_smc %>%
    select(-est_median:-est_ltail) %>%
    pivot_longer(rhat_median:rmse_ltail, names_to=c("var", "stat"), names_sep="_") %>%
    mutate(value = case_when( # variable transformations
        var == "sd" & value == 0 ~ NA_real_,
        var == "rhat" ~ value - 1,
        TRUE ~ value
    ))
# write_rds(res, here("validation/out/6x6/gsmc_valid_6x6in6.rds"), compress="gz")
# write_rds(res, here("~/Desktop/gsmc_valid_6x6in6.rds"), compress="gz")
# quit()

## diagnostic plots -------

# library(here)
# devtools::load_all(recompile=F)
# library(sf)
# library(tidyverse)
# library(patchwork)
# library(scales)
# res = read_rds(here("~/Desktop/gsmc_valid_6x6in6.rds"))

make_iter_plot = function(x, ylab=NULL, ymax=NA, ymin=NA) {
    filter(res, var==x) %>%
        mutate(stat = fct_inorder(
            c(median="Pr(edges rem. = 22)",
              htail="Pr(edges rem. > 26)",
              ltail="Pr(edges rem. = 18)")[stat])
        ) %>%
    ggplot(aes(iter, value,  lty=stat)) +
        geom_hline(yintercept=0.0, color="#00000077") +
        geom_line(linewidth=0.5, color="#444444") +
        geom_point(size=1.5, color="#444444") +
        scale_x_continuous(trans="log10", labels=comma) +
        labs(x="Sample size", y=NULL, title=ylab, lty="Statistic",
             color="Algorithm", shape="Algorithm") +
        coord_cartesian(ylim=c(ymin, ymax)) +
        { if (x == "rhat") scale_y_continuous(labels=function(y) number(1+y, 0.01)) } +
        # { if (x %in% c("sd", "rmse")) scale_y_continuous(trans="sqrt") } +
        theme_bw(base_family="Times", base_size=10) +
        theme(plot.margin=margin(r=10),
              legend.margin=margin(),
              plot.background=element_blank(),
              panel.background=element_blank())
}

p = make_iter_plot("rhat", expression(hat(R))) +
    make_iter_plot("sd", "Standard error", ymax=0.1) +
    make_iter_plot("bias", "Bias", ymax=0.164, ymin=-0.164) +
    make_iter_plot("rmse", "RMSE", ymax=0.1) +
plot_layout(nrow=1, guides="collect") &
    theme(legend.position="bottom",
          plot.margin=unit(rep(1, 4), "mm"))

ggsave(here("~/Desktop//gsmc_diagn_6x6.pdf"), plot=p, width=8, height=3.25)

# Check validity ---------
make_valid_plot = function(M, conf=0.9, by_chain=TRUE) {
    pl_smc = pls_smc[[as.character(M)]]

    # edge cut compactness
    edge_rg = do.call(seq, as.list(range(pl_enum$edges)))
    d_hist = map_dfr(edge_rg, function(k) {
        redist_ci(pl_smc, edges == k, 1, conf=conf, by_chain=by_chain) %>%
            suppressWarnings() %>%
            `colnames<-`(c("est", "low", "high")) %>%
            mutate(edges = k,
                   true = weighted.mean(by_plan(pl_enum$edges, n) == k, w))
    })

    ggplot(d_hist, aes(edges)) +
        geom_col(aes(y=true), fill=PAL[1]) +
        geom_pointrange(aes(y=est, ymin=low, ymax=high), color="#444444",
                        linewidth=0.7, fatten=1.5) +
        scale_y_continuous(str_glue("Proportion of plans"), labels=percent,
                           expand=expansion(mult=c(0, 0.04))) +
        coord_cartesian(ylim=c(0, NA)) +
        labs(title=str_glue("{comma(M)} samples per run"), x="Number of edges removed") +
        theme_bw(base_family="Times", base_size=10)
}


p = make_valid_plot(100) + make_valid_plot(1000) + make_valid_plot(10000) +
    plot_layout(nrow=1) & theme(plot.margin=margin(l=0.8, unit="mm"))
ggsave(here("~/Desktop/gsmc_valid_6x6.pdf"), plot=p, width=8, height=2.5)
