library(tidyverse)

data(iowa)
ia <- redist_map(iowa, existing_plan=cd_2010, pop_tol=0.01)
plans <- redist_smc(ia, 1000, adapt_k_thresh=1, adjust_labels=TRUE)

get_sampling_info(plans)$diagnostics$est_k
devs = map_dfr(get_sampling_info(plans)$diagnostics$adapt_devs, function(x) {
    tibble(k=seq_along(x), n_ok=x)
}, .id="iter") %>%
    mutate(iter=as.integer(iter))

filter(devs, k<10) %>%
ggplot(aes(k, log(n_ok), color=iter, group=iter)) +
    geom_line(size=1.0)

filter(devs, iter==2, draw==8) %>%
    qplot(dev, data=.)
