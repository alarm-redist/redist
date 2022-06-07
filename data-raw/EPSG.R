library(tidyverse)

# reproducible code for making EPSG lookup
make_epsg_table <- function() {
    raw <- as_tibble(rgdal::make_EPSG()) %>%
        select(code, note)
    state_regex <- paste0("(", paste0(datasets::state.name, collapse = "|"), ")")
    epsg_regex <- str_glue("NAD83(\\(HARN\\))? / {state_regex} ?[A-Za-z. ]*$")
    epsg_d <- filter(raw, (code > 2500L & code < 2900L) | (code > 3300L & code < 3400L),
        str_detect(note, epsg_regex)) %>%
        mutate(state = str_match(note, epsg_regex)[, 3],
            priority = str_detect(note, "HARN") + str_detect(note, "Central")) %>%
        group_by(state) %>%
        arrange(desc(priority)) %>%
        slice(1) %>%
        ungroup() %>%
        select(code, state) %>%
        rows_insert(tibble(code = 2784L, state = "Hawaii"), by = "state") %>%
        arrange(state)

    codes <- as.list(epsg_d$code)
    names(codes) <- datasets::state.abb
    codes
}

EPSG <- make_epsg_table()

usethis::use_data(EPSG, overwrite = TRUE)
