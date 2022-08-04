library(tidyverse)
library(callr)
library(microbenchmark)

devtools::load_all()

data(iowa)
ia <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01)

x = redist_smc(ia, 100, ncores=1L, silent=TRUE)
redist.plot.plans(x, 1:4, ia)

microbenchmark(redist_smc(ia, 500, ncores=1L, silent=TRUE), times=10)
# 100 -> 182ms
# 500 -> 976ms

va = readRDS("../partisan-bias/data-raw/virginia/VA_map.rds") |>
    set_pop_tol(0.01)

x = redist_smc(va, 500, ncores=1L, verbose=TRUE)



########################

res = callr::r_bg(function(...) { redist::redist_smc(...) },
                args=list(map=ia, nsims=1000, ncores=1L, silent=TRUE),
                poll_connection=FALSE)

d = tibble(raw = system2("sample", args=as.character(res$get_pid()), stdout=TRUE))

write_lines(d$raw, "~/Desktop/sample.txt")
d = tibble(raw = read_lines("~/Desktop/sample.txt"))

d |>
    filter(str_starts(str_trim(raw), "\\+")) |>
    mutate(n_mark = str_count(raw, "([+|!]|: )"),
           raw = str_remove_all(raw, "([+|!]|: )"),
           raw = str_replace_all(raw, "[()]", " "),
           raw = str_squish(raw)) |>
    filter(str_detect(raw, "redist\\.so")) |>
    mutate(time = as.numeric(word(raw, 1)),
           fn = word(raw, 2),
           loc = word(raw, -1)) |>
    # filter(!str_starts(loc, "(main|eval)\\.c")) |>
    filter(!str_starts(fn, "std::")) |>
    select(-raw) |>
    group_by(fn) |>
    summarize(time = sum(time),
              depth = median(n_mark)) |>
    mutate(time = time / max(time)) |>
    arrange(desc(time)) |>
    as.data.frame()
