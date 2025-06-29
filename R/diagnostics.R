
#' Diagnostic information on sampled plans
#'
#' Prints diagnostic information, which varies by algorithm. All algorithms
#' compute the [plans_diversity()] of the samples.
#'
#' For SMC and MCMC, if there are multiple runs/chains, R-hat values will be
#' computed for each summary statistic. These values should be close to 1.
#' If they are not, then there is too much between-chain variation, indicating
#' that there are not enough samples. R-hat values are calculated after
#' rank-normalization and folding.  MCMC chains are split in half before R-hat
#' is computed. For summary statistics that vary across districts, R-hat is
#' calculated for the first district only.
#'
#' For SMC, diagnostics statistics include:
#'
#' * **Effective samples**: the effective sample size at each iteration, computed
#' using the SMC weights. Larger is better. The percentage in parentheses is the
#' ratio of the effective samples to the total samples.
#' * **Acceptance rate**: the fraction of drawn spanning trees which yield a valid
#' redistricting plan within the population tolerance. Very small values (< 1%)
#' can indicate a bottleneck and may lead to a lack of diversity.
#' * **Standard deviation of the log weights**: More variable weights (larger s.d.)
#' indicate less efficient sampling. Values greater than 3 are likely problematic.
#' * **Maximum unique plans:** an upper bound on the number of unique redistricting
#' plans that survive each stage. The percentage in parentheses is the ratio of
#' this number to the total number of samples. Small values (< 100) indicate a
#' bottleneck, which leads to a loss of sample diversity and a higher variance.
#' * **Estimated `k` parameter**: How many spanning tree edges were considered for
#' cutting at each split. Mostly informational, though large jumps may indicate
#' a need to increase `adapt_k_thresh`.
#' * **Bottleneck**: An asterisk will appear in the right column if a bottleneck
#' appears likely, based on the values of the other statistics.
#'
#' In the event of problematic diagnostics, the function will provide
#' suggestions for improvement.
#'
#' @param object a [redist_plans] object
#' @param district For R-hat values, which district to use for district-level
#' summary statistics. We strongly recommend calling `match_numbers()` or
#' `number_by()` before examining these district-level statistics. A value of
#' `FALSE` indicates to compute the rhats for all districts
#' @param all_runs When there are multiple SMC runs, show detailed summary
#' statistics for all runs (the default), or only the first run?
#' @param vi_max The maximum number of plans to sample in computing the pairwise
#' variation of information distance (sample diversity).
#' @param use_order_stats Whether or not to compute rhats on the ordered district
#' statistics.
#' @param \dots additional arguments (ignored)
#'
#' @return A data frame containing diagnostic information, invisibly.
#'
#' @examples
#' data(iowa)
#' iowa_map <- redist_map(iowa, ndists = 4, pop_tol = 0.1)
#' plans <- redist_smc(iowa_map, 100)
#' summary(plans)
#'
#' @method summary redist_plans
#' @concept analyze
#' @md
#' @export
summary.redist_plans <- function(
        object, district = FALSE, all_runs = TRUE, vi_max = 100,
        use_order_stats = TRUE,...) {
    cli::cli_process_done(done_class = "") # in case an earlier

    algo <- attr(object, "algorithm")
    name <- deparse(substitute(object))

    object <- subset_sampled(object)
    all_diagn <- attr(object, "diagnostics")
    all_run_info <- attr(object, "run_information")
    plans_m <- get_plans_matrix(object)
    n_samp <- ncol(plans_m)
    n_distr <- attr(object, "ndists")
    if (is.null(n_distr)) n_distr <- max(plans_m[, 1])

    fmt_comma <- function(x) format(x, nsmall = 0, digits = 1, big.mark = ",")
    if (n_distr == 1 || nrow(plans_m) == 1) {
        cli_text("{fmt_comma(n_samp)}{cli::qty(n_samp)} sampled plan{?s} of
                 {n_distr} district{?s} on
                 {fmt_comma(nrow(plans_m))}{cli::qty(nrow(plans_m))} unit{?s}")
        return(invisible(1))
    }

    prec_pop <- attr(object, "prec_pop")
    if (is.null(prec_pop)) {
        cli_warn(c("Precinct population missing; plan diversity estimates may be misleading.",
                   ">" = 'Run `attr({name}, "prec_pop") <- <map object>$<pop column>` to fix.'))
        prec_pop <- rep(1, nrow(plans_m))
    }
    est_div <- plans_diversity(object, total_pop = prec_pop, n_max = vi_max)
    div_rg <- format(quantile(est_div, c(0.1, 0.9)), digits = 2)
    div_bad <- (mean(est_div) <= 0.35) || (mean(est_div <= 0.25) > 0.1)

    summary_supported_algs <- c(SMC_ALG_TYPE, MS_SMC_ALG_TYPE, MCMC_ALG_TYPE, "flip")

    # ignore if not a supported algorithm
    if(!algo %in% summary_supported_algs){
        cli_abort("{.fn summary} is not supported for the {toupper(algo)} algorithm.")
        return(invisible(1))
    }

    # bool for not flip
    revamped_alg <- algo %in% c(SMC_ALG_TYPE, MS_SMC_ALG_TYPE, MCMC_ALG_TYPE)

    # get the display name
    alg_display_name <- case_when(
        algo == SMC_ALG_TYPE ~ "SMC",
        algo == MS_SMC_ALG_TYPE ~ "SMC with Merge-Split MCMC Steps",
        algo == MCMC_ALG_TYPE ~ "Merge-Split MCMC",
        algo == "flip" ~ "Flip MCMC"
    )
    cli_text("{.strong {alg_display_name}:} {fmt_comma(n_samp)} sampled plans of {n_distr}
                 districts on {fmt_comma(nrow(plans_m))} units")

    # if revamped alg check that sampling space and splitting methods all the same
    if(revamped_alg){
        # check same sampling space
        all_sampling_spaces <- sapply(all_run_info, function(x) x$sampling_space)
        if(length(unique(all_sampling_spaces)) != 1){
            cli::cli_abort("{.fn summary} is not supported for plans sampled using different sampling spaces.")
            return(invisible(1))
        }
        sampling_space <- all_sampling_spaces[1]
        display_sampling_space <- case_when(
            sampling_space == GRAPH_PLAN_SPACE_SAMPLING ~ "Graph Space",
            sampling_space == FOREST_SPACE_SAMPLING ~ "Spanning Forest Space",
            sampling_space == LINKING_EDGE_SPACE_SAMPLING ~ "Linking Edge Forest Space"
        )
        # check same splitting method
        all_splitting_methods <- sapply(all_run_info, function(x) x$split_method)
        if(length(unique(all_splitting_methods)) != 1){
            cli::cli_abort("{.fn summary} is not supported for plans sampled using different splitting methods")
            return(invisible(1))
        }
        split_method <- all_splitting_methods[1]
        display_splitting_method <- case_when(
            split_method == NAIVE_K_SPLITTING ~ "Naive Top K",
            split_method == UNIF_VALID_EDGE_SPLITTING ~ "Uniform Valid Edge",
            split_method == EXP_BIGGER_ABS_DEV_SPLITTING ~ "Exponential Absolute Deviance"
        )

        # check the splitting parameters are all the same
        all_split_params <- lapply(all_diagn, function(x) x$split_params)
        for (i in seq_len(length(all_split_params))) {
            if(!identical(all_split_params[[1]], all_split_params[[i]])){
                cli::cli_abort("{.fn summary} is not supported for plans sampled using different splitting parameters")
                return(invisible(1))
            }
        }
        split_params <- all_split_params[[1]]

        cli::cli_text("Plans sampled on {display_sampling_space} using the {display_splitting_method} forward kernel.")

        if(split_method == NAIVE_K_SPLITTING){
            # only display adapt k threshold if k values estimated
            if(split_params$estimate_cut_k){
                cli::cli_text("Forward kernel parameters: {.arg adapt_k_thresh}={format(split_params$adapt_k_thresh, digits=3)}")
            }
        }else if(split_method == UNIF_VALID_EDGE_SPLITTING){

        }else if(split_method == EXP_BIGGER_ABS_DEV_SPLITTING){
            cli::cli_text("Forward Kernel Parameters: {.arg splitting_alpha}={format(split_params$splitting_alpha, digits=2)}")
        }
    }


    # print algorithm specific parameters
    if (algo %in% c(SMC_ALG_TYPE, MS_SMC_ALG_TYPE)) {
        cli::cli_text("SMC Parameters: {.arg weight_type}={format(all_run_info[[1]]$weight_type, digits=2)}
        \u2022 {.arg seq_alpha}={format(all_diagn[[1]]$seq_alpha, digits=2)}
        \u2022 {.arg pop_temper}={format(all_diagn[[1]]$pop_temper, digits=3)}")
    }else if(algo == MCMC_ALG_TYPE){
        cli::cli_text(
        "MCMC Parameters: {.arg warmup}={format(all_diagn[[1]]$warmup)}
        \u2022 {.arg thin}={format(all_diagn[[1]]$thin)}
        \u2022 {.arg thin}={format(all_diagn[[1]]$total_steps)}")
    }


    cli_text("Plan diversity 80% range: {div_rg[1]} to {div_rg[2]}")
    if (div_bad) cli::cli_alert_danger("{.strong WARNING:} Low plan diversity")


    # now compute rhats if more than 1 chain
    one_district_only <- 1 <= district && district <= n_distr
    cols <- names(object)
    addl_cols <- setdiff(cols, c("chain", "draw", "district", "total_pop", "nseats", "mcmc_accept"))
    if(one_district_only){
        idx <- seq_len(n_samp)
        if ("district" %in% cols) idx <- as.integer(district) + (idx - 1)*n_distr

        const_cols <- vapply(addl_cols, function(col) {
            x <- object[[col]][idx]
            all(is.na(x)) || all(x == x[1]) ||
                any(tapply(x, object[['chain']][idx], FUN = function(z) length(unique(z))) == 1)
        }, numeric(1))
    }else{
        const_cols <- vapply(addl_cols, function(col) {
            x <- object[[col]]
            all(is.na(x)) || all(x == x[1]) ||
                any(tapply(x, object[['chain']], FUN = function(z) length(unique(z))) == 1)
        }, numeric(1))
    }
    addl_cols <- addl_cols[!const_cols]

    # do nothing if no additional columns or no chain column
    if(length(addl_cols) > 0 && "chain" %in% cols){
        rhats_computed <- TRUE
        split_rhat = algo %in% c(MCMC_ALG_TYPE, "flip")

        if(1 <= district && district <= n_distr ){
            rhat_df <- object |>
                filter(!is.na(chain) & district == !!district) |>
                pivot_longer(
                    cols = all_of(addl_cols),
                    names_to = "variable",
                    values_to = "value"
                ) |>
                filter(!is.na(value)) |>
                select(variable, chain, district, value) |>
                group_by(variable, district) |>
                summarise(rhat = diag_rhat(value, chain, split=split_rhat), .groups = "drop")
        }
        if(use_order_stats){
            # else compute the rhats on the order statistics
            rhat_df <- object |>
                filter(!is.na(chain)) |>
                pivot_longer(
                    cols = all_of(addl_cols),
                    names_to = "variable",
                    values_to = "value"
                ) |>
                filter(!is.na(value)) |>
                group_by(variable, chain, draw) |>
                mutate(district = row_number(value)) |>
                ungroup() |>
                select(variable, chain, district, value) |>
                group_by(variable, district) |>
                summarise(rhat = diag_rhat(value, chain, split=split_rhat), .groups = "drop")
        }else{
            # else compute rhats using existing district number
            rhat_df <- object |>
                filter(!is.na(chain)) |>
                pivot_longer(
                    cols = all_of(addl_cols),
                    names_to = "variable",
                    values_to = "value"
                ) |>
                filter(!is.na(value)) |>
                select(variable, chain, district, value) |>
                group_by(variable, district) |>
                summarise(rhat = diag_rhat(value, chain, split=split_rhat), .groups = "drop")
        }
        ordered_str <- ifelse(use_order_stats, "ordered ", "")
        cli::cli_text("Largest R-hat values for {ordered_str}summary statistics:\n")
        max_rhats_df <- rhat_df |>
            group_by(variable) |>
            summarise(max_rhat = max(rhat))
        max_rhats <- max_rhats_df |>
            pull(max_rhat)
        names(max_rhats) <- max_rhats_df |>
            pull(variable)
        rhats_p <- vapply(max_rhats, function(x){
            ifelse(x < 1.05, sprintf('%.3f', x), paste0('\U274C', round(x, 3)))
        }, FUN.VALUE = character(1))
        print(noquote(rhats_p))

        if (any(rhat_df$rhat >= 1.05)){
            warn_converge <- TRUE
            cli::cli_alert_danger("{.strong WARNING:} Chains have not converged.")
        }
    }else{
        warn_converge <- FALSE
        rhats_computed <- FALSE
    }




    # Now print algorithm specific diagnostics
    if (algo %in% c(SMC_ALG_TYPE, MS_SMC_ALG_TYPE)) {

        run_dfs <- list()
        n_runs <- length(all_diagn)
        warn_bottlenecks <- FALSE

        for (i in seq_len(n_runs)) {
            diagn <- all_diagn[[i]]
            n_samp <- nrow(diagn$ancestors)



            # modified to deal with the fact I changed plan data at some point

            # if new do nothing
            if(length(diagn$sd_lp) == length(diagn$nunique_original_ancestors)){
                a_unique_orig_ancestors <- diagn$nunique_original_ancestors
            }else{
                a_unique_orig_ancestors <- c(diagn$nunique_original_ancestors, NA)
            }

            run_dfs[[i]] <- tibble(n_eff = c(diagn$step_n_eff, diagn$n_eff),
                                   eff = c(diagn$step_n_eff, diagn$n_eff)/n_samp,
                                   accept_rate = c(diagn$accept_rate, NA),
                                   sd_log_wgt = diagn$sd_lp,
                                   max_unique = diagn$unique_survive,
                                   est_k = c(diagn$est_k, NA),
                                   unique_original = a_unique_orig_ancestors)

            tbl_print <- as.data.frame(run_dfs[[i]])
            min_n <- max(0.05*n_samp, min(0.4*n_samp, 100))
            bottlenecks <- dplyr::coalesce(with(tbl_print, pmin(max_unique, n_eff) < min_n), FALSE)
            warn_bottlenecks <- warn_bottlenecks || any(bottlenecks)
            tbl_print$bottleneck <- ifelse(bottlenecks, " * ", "")
            tbl_print$n_eff <- with(tbl_print,
                                    str_glue("{fmt_comma(n_eff)} ({sprintf('%0.1f%%', 100*eff)})"))
            tbl_print$eff <- NULL
            tbl_print$accept_rate <- with(tbl_print, sprintf("%0.1f%%", 100*accept_rate))
            max_pct <- with(tbl_print, max_unique/(-n_samp * expm1(-1)))
            tbl_print$max_unique <- with(tbl_print,
                                         str_glue("{fmt_comma(max_unique)} ({sprintf('%3.0f%%', 100*max_pct)})"))

            names(tbl_print) <- c("Eff. samples (%)", "Acc. rate",
                                  "Log wgt. sd", " Max. unique",
                                  "Est. k", "Unique Original Ancestors" , "")
            rownames(tbl_print) <- c(paste("Split", seq_len(nrow(tbl_print) - 1)), "Resample")

            if (i == 1 || isTRUE(all_runs)) {
                cli_text("Sampling diagnostics for SMC run {i} of {n_runs} ({fmt_comma(n_samp)} samples)")
                print(tbl_print, digits = 2)
                cat("\n")
            }
        }
        out <- bind_rows(run_dfs)

        cli::cli_li(cli::col_grey("
            Watch out for low effective samples, very low acceptance rates (less than 1%),
            large std. devs. of the log weights (more than 3 or so),
            and low numbers of unique plans.
            R-hat values for summary statistics should be between 1 and 1.05."))

        if (div_bad) {
            cli::cli_li("{.strong Low diversity:} Check for potential bottlenecks.
                        Increase the number of samples.
                        Examine the diversity plot with
                        `hist(plans_diversity({name}), breaks=24)`.
                        Consider weakening or removing constraints, or increasing
                        the population tolerance. If the acceptance rate drops
                        quickly in the final splits, try increasing
                        {.arg pop_temper} by 0.01.")
        }
        if (warn_converge) {
            cli::cli_li("{.strong SMC convergence:} Increase the number of samples.
                        If you are experiencing low plan diversity or bottlenecks as well,
                        address those issues first.")
        }
        if (warn_bottlenecks) {
            cli::cli_li("(*) {.strong Bottlenecks found:} Consider weakening or removing
                        constraints, or increasing the population tolerance.
                        If the acceptance rate drops quickly in the final splits,
                        try increasing {.arg pop_temper} by 0.01.
                        If the weight variance (Log wgt. sd) increases steadily
                        or is particularly large for the \"Resample\" step,
                        consider increasing {.arg seq_alpha}.
                        To visualize what geographic areas may be causing problems,
                        try running the following code. Highlighted areas are
                        those that may be causing the bottleneck.\n\n")
            code <- str_glue("plot(<map object>, rowMeans(as.matrix({name}) == <bottleneck iteration>))")
            cli::cat_line("    ", cli::code_highlight(code, "Material"))
        }
    }else if (algo %in% c("mergesplit", 'flip')) {

        accept_rate <- sprintf("%0.1f%%", 100*attr(object, "mh_acceptance"))
        cli_text("Chain acceptance rate{?s}: {accept_rate}")

        if(rhats_computed){
            out <- max_rhats_df
        }else{
            out <- tibble(accept_rate = attr(object, "mh_acceptance"),
                          div_q10 = div_rg[1],
                          div_q90 = div_rg[2])
        }

        cli::cli_li(cli::col_grey("
            Watch out for low acceptance rates (less than 10%).
            R-hat values for summary statistics should be between 1 and 1.05.
            For district-level statistics (like district partisan leans), you
            should call `match_numbers()` or `number_by()` before examining
            the R-hat values."))

        if (div_bad) {
            cli::cli_li("{.strong Low diversity:} Increase the number of samples.
                        Examine the diversity plot with
                        `hist(plans_diversity({name}), breaks=24)`.
                        Consider weakening or removing constraints, or increasing
                        the population tolerance.")
        }
        if (warn_converge) {
            cli::cli_li("{.strong Chain convergence:} Increase the number of samples.
                        If you are experiencing low plan diversity, address that issue first.")
        }
    } else {
        cli_abort("{.fn summary} is not supported for the {toupper(algo)} algorithm.")
    }

    invisible(out)
}


diag_fold <- function(x) {
    abs(x - median(x))
}

diag_ranknorm <- function(x) {
    qnorm(rank(x)/(length(x) + 1))
}

diag_calc_rhat <- function(x, grp) {
    n <- mean(table(grp))
    var_between <- n*var(tapply(x, grp, mean))
    var_within <- mean(tapply(x, grp, var))
    sqrt((var_between/var_within + n - 1)/n)
}

diag_rhat <- function(x, grp, split = FALSE) {
    if (split) {
        lengths <- rle(grp)$lengths
        grp <- grp + do.call(c, lapply(lengths, function(l) rep(c(0.0, 0.5), each = l/2)))
    }

    max(diag_calc_rhat(diag_ranknorm(x), grp),
        diag_calc_rhat(diag_ranknorm(diag_fold(x)), grp))
}



#' Get k-step Ancestors of particles
#'
#' Tells you for a given index what its original ancestor was.
#'
#' @param parent_mat Ancestor matrix where entry [i,j] equals the index of the
#' parent of particle i after step j
#' @param steps_back How many steps back we should find the ancestors of (so
#' parents are 1 step ancestors). Defaults to going all the way back to the
#' original ancestors
#' @param start_col Which column to start at. Defaults to the last column
#'
#' @returns indices of k-step ancestors for particles at iteration `start_col`
#' @md
#' @export
get_k_step_ancestors <- function(parent_mat, steps_back = NULL, start_col = NULL){
    # check the matrix is not zero indexed
    assertthat::assert_that(
        ! 0 %in% parent_mat,
        msg = "parent_mat must be 1-indexed, not 0 indexed!"
    )

    nparent_cols <- ncol(parent_mat)

    # if null then just start on final set of plans
    if(is.null(start_col)){
        start_col <- nparent_cols
    }else{
        # else assert starting column is between 1 and number of cols
        assertthat::assert_that(
            assertthat::is.scalar(start_col) &&
                start_col <= nparent_cols &&
                start_col > 1,
            msg = sprintf("Error!\nInput start_col=`%s` is not valid.
                          Input must be a number between 1 and %d",
                          toString(start_col), nparent_cols)
        )
    }

    # check steps back is between [1, nnparent_cols -1]
    if(is.null(steps_back)){
        steps_back <- nparent_cols-1
    }else{
        assertthat::assert_that(
            assertthat::is.scalar(steps_back) &&
                steps_back <= start_col-1 &&
                steps_back >= 1,
            msg = sprintf(
                "Error!\nInput steps_back=`%s` is not valid.
Input must be between 1 and the start_col value (you input %d)",
                toString(steps_back), steps_back)
        )
    }

    # vector where index i maps to the index of its ancestor
    # initialize to this
    ancestor <- 1:nrow(parent_mat)


    # iterate through each step back we select the successive parent indices
    for (j in start_col:(start_col - steps_back + 1)) {
        ancestor <- parent_mat[ancestor, j]
    }
    return(ancestor)
}



#' Get Original Ancestor Matrix of particles
#'
#' Gets a matrix of the original ancestors (ie first splits) of the particles
#' at each steo,
#'
#' @param parent_mat Ancestor matrix where entry [i,j] equals the index of the
#' parent of particle i after step j
#'
#' @returns indices of originals ancestors for particles at every iteration. So
#' entry [s,j] is the original ancestor of particle j after step s
#' @md
#' @export
get_original_ancestors_mat <- function(parent_mat){

    #if the matrix has one column then its just itself
    if(ncol(parent_mat) == 1){
        return(parent_mat)
    }

    # get the original ancestors at every step from the parent matrix
    original_ancestor_mat <- sapply(
        2:ncol(parent_mat),
        function(col_num) redist::get_k_step_ancestors(parent_mat, steps_back = col_num-1, start_col = col_num)
    )

    # add the first column where every particles ancestor is iteslf
    original_ancestor_mat <- cbind(
        1:nrow(parent_mat),
        original_ancestor_mat
    )

    return(original_ancestor_mat)
}
