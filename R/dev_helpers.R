# THIS FILE IS NOT INCLUDED IN BUILDS

# call ld_ia() to create a redist_map and redist_plans object for Iowa
ld_ia <- function() {
    data(iowa)
    ia <<- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01)
    plans <<- redist_smc(ia, 100, silent = TRUE)
}

enforce_style <- function() {
    R_style <- function(...) {
        x <- styler::tidyverse_style(scope = I(c("spaces", "indention", "tokens")),
            indent_by = 4,
            strict = FALSE,
            start_comments_with_one_space = TRUE,
            math_token_spacing = styler::specify_math_token_spacing(
                zero = c("'^'", "'*'", "'/'"),
                one = c("'+'", "'-'")))
        x
    }

    styler::cache_activate()
    styler::style_pkg(style = R_style,
                      exclude_files=c("R/redist_smc.R", "R/redist_ms.R", "R/redist_ms_parallel.R"),
                      exclude_dirs = "explore")
}
