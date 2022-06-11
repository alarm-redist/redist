#' 'Flip' Markov Chain Monte Carlo Redistricting Simulation
#'
#' \code{redist_flip} provides a tidy interface to the methods in
#' \code{\link{redist.flip}}.
#'
#' This function allows users to simulate redistricting plans
#' using the Markov Chain Monte Carlo methods of Fifield et al. Several
#' constraints corresponding to substantive requirements in the redistricting
#' process are implemented, including population parity and geographic
#' compactness. In addition, the function includes multiple-swap and simulated
#' tempering functionality to improve the mixing of the Markov Chain.
#'
#' \code{redist_flip} allows for Gibbs constraints to be supplied via a list object
#' passed to \code{constraints}. This is a change from the original \code{redist.flip}
#' behavior to allow for a more straightforward function call when used within a pipe.
#' A key difference between \code{redist_flip} and \code{redist.flip} is that
#' \code{redist_flip} uses a small compactness constraint by default, as this improves
#' the realism of the maps greatly and also leads to large speed improvements.
#' (One of the most time consuming aspects of the flip MCMC backend is checking for
#' district shattering, which is slowed down even further by non-compact districts.
#' As such, it is recommended that all flip simulations use at least a minimal compactness
#' constraint, even if you weaken it from the default settings.) The default is
#' a \code{compact} constraint using the \code{edges-removed} metric with a
#' weight of 0.6. For very small maps (< 100 precincts), you will likely want to
#' weaken (lower) this constraint, while for very large maps (> 5000 precincts),
#' you will likely want to strengthen (increase) this constraint. Otherwise,
#' for most maps, the default constraint should be a good starting place.
#'
#' \code{redist_flip} samples from a known target distribution which can be described
#' using the \code{constraints}. The following describes the constraints available. The general
#' advice is to set weights in a way that gets between 20\% and 40\% acceptance
#' on average, though more tuning advice is available in the vignette on using
#' MCMC methods.Having too small of an acceptance rate indicates that the weights
#' within \code{constraints} are too large and will impact sampling efficiency.
#' If the Metropolis Hastings acceptance rate is too large, this may impact the
#' target distribution, but may be fine for general exploration of possible maps.
#'
#' There are currently 9 implemented constraint types, though `\code{compact} and
#' \code{partisan} have sub-types which are specified via a character \code{metric}
#' within their respective list objects. The constraints are as follows:
#' * \code{compact} - biases the algorithm towards drawing more compact districts.
#' * weight - the coefficient to put on the Gibbs constraint
#' * metric - which metric to use. Must be one of \code{edges-removed} (the default),
#' \code{polsby-popper}, \code{fryer-holden}, or \code{log-st}. Using Polsby Popper
#' is generally not recommended, as \code{edges-removed} is faster and highly correlated.
#' \code{log-st} can be used to match the target distribution of \code{redist_smc} or
#' \code{redist_mergesplit}.
#' * areas - Only used with \code{polsby-popper} - A vector of precinct areas.
#' * borderlength_mat - Only used with \code{polsby-popper} - A matrix of precinct
#' border lengths.
#' * ssdmat - Only used with \code{fryer-holden} - A matrix of squared distances between
#' precinct centroids.
#' * ssd_denom - Only used with \code{fryer-holden} - a positive integer to use
#' as the normalizing constant for the Relative Proximity Index.
#' * \code{population} - A Gibbs constraint to complement the hard population
#' constraint set by \code{pop_tol}. This penalizes moves which move away from smaller
#' population parity deviations. It is very useful when an \code{init_plan} sits
#' outside of the desired \code{pop_tol} but there are substantive reasons to use
#' that plan. This constraint uses the input to \code{total_pop}.
#' * weight - the coefficient to put on the Gibbs constraint
#' * \code{countysplit} This is a Gibbs constraint to minimize county splits. Unlike
#' SMC's county constraint, this allows for more than \code{ndists - 1} splits and
#' does not require that counties are contiguous.
#' * weight - the coefficient to put on the Gibbs constraint
#' * \code{hinge} This uses the proportion of a group in a district and matches to the
#' nearest target proportion, and then creates a penalty of
#' \eqn{\sqrt{max(0, nearest.target - group.pct)}}.
#' * weight - the coefficient to put on the Gibbs constraint
#' * minorityprop - A numeric vector of minority proportions (between 0 and 1) which
#' districts should aim to have
#' * \code{vra} This takes two target proportions of the presence of a minority group
#' within a district. \eqn{(|target.min - group.pct||target.other - group.pct|)^{1.5})}
#' * weight - the coefficient to put on the Gibbs constraint
#' * target_min - the target minority percentage. Often, this is set to 0.55 to encourage
#' minority majority districts.
#' * target_other - the target minority percentage for non majority minority districts.
#' * \code{minority} This constraint sorts the districts by the proportion of a group in
#' a district and compares the highest districts to the entries of minorityprop.
#' This takes the form \eqn{\sum_{i=1}^{n} \sqrt{|group.pct(i) - minorityprop(i)| }} where n
#' is the length of minorityprop input.
#' * weight - the coefficient to put on the Gibbs constraint
#' * minorityprop - A numeric vector of minority proportions (between 0 and 1) which
#' districts should aim to have
#' * \code{similarity} This is a status-quo constraint which penalizes plans which
#' are very different from the starting place. It is useful for local exploration.
#' * weight - the coefficient to put on the Gibbs constraint
#' * \code{partisan} This is a constraint which minimizes partisan bias, either as
#' measured as the difference from proportional representation or as the magnitude of
#' the efficiency gap.
#' * weight - the coefficient to put on the Gibbs constraint
#' * rvote - An integer vector of votes for Republicans or other party
#' * dvote - An integer vector of votes for Democrats or other party
#' * metric - which metric to use. Must be one of \code{proportional-representation}
#' or \code{efficiency-gap}.
#' * \code{segregation} This constraint attempts to minimize the degree of dissimilarity
#' between districts by group population.
#' * weight - the coefficient to put on the Gibbs constraint
#'
#'
#' @param map A \code{\link{redist_map}} object.
#' @param nsims The number of samples to draw, not including warmup.
#' @param warmup The number of warmup samples to discard.
#' @param init_plan A vector containing the congressional district labels
#' of each geographic unit. The default is \code{NULL}. If not provided,
#' a random initial plan will be generated using \code{redist_smc}. You can also
#' request to initialize using \code{redist.rsg} by supplying 'rsg', though this is
#' not recommended behavior.
#' @param constraints A `redist_constr` object.
#' @param nthin The amount by which to thin the Markov Chain. The
#' default is \code{1}.
#' @param eprob The probability of keeping an edge connected. The
#' default is \code{0.05}.
#' @param lambda lambda The parameter determining the number of swaps to attempt
#' each iteration of the algorithm. The number of swaps each iteration is
#' equal to Pois(\code{lambda}) + 1. The default is \code{0}.
#' @param temper Whether to use simulated tempering algorithm. Default is FALSE.
#' @param betaseq Sequence of beta values for tempering. The default is
#' \code{powerlaw} (see Fifield et. al (2020) for details).
#' @param betaseqlength Length of beta sequence desired for
#' tempering. The default is \code{10}.
#' @param betaweights betaweights Sequence of weights for different values of
#' beta. Allows the user to upweight certain values of beta over
#' others. The default is \code{NULL} (equal weighting).
#' @param adapt_lambda adapt_lambda Whether to adaptively tune the lambda parameter so that the Metropolis-Hastings
#' acceptance probability falls between 20% and 40%. Default is FALSE.
#' @param adapt_eprob eprob Whether to adaptively tune the edgecut probability parameter so that the
#' Metropolis-Hastings acceptance probability falls between 20% and 40%. Default is
#' FALSE.
#' @param exact_mh  Whether to use the approximate (FALSE) or exact (TRUE)
#' Metropolis-Hastings ratio calculation for accept-reject rule. Default is FALSE.
#' @param adjswaps Flag to restrict swaps of beta so that only
#' values adjacent to current constraint are proposed. The default is
#' \code{TRUE}.
#' @param init_name a name for the initial plan, or \code{FALSE} to not include
#' the initial plan in the output.  Defaults to the column name of the
#' existing plan, or "\code{<init>}" if the initial plan is sampled.
#' @param verbose Whether to print initialization statement. Default is \code{TRUE}.
#'
#' @return A \code{\link{redist_plans}} object containing the simulated plans.
#' @concept simulate
#' @export
#'
#' @references Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander
#' Tarr. (2016) "A New Automated Redistricting Simulator Using Markov Chain Monte
#' Carlo." Working Paper. Available at
#' \url{http://imai.princeton.edu/research/files/redist.pdf}.
#'
#' @importFrom dplyr bind_cols
#' @importFrom rlang eval_tidy enquo
#' @importFrom utils capture.output
#' @md
#' @examples
#' data(iowa)
#' iowa_map <- redist_map(iowa, ndists = 4, existing_plan = cd_2010, total_pop = pop,
#'     pop_tol = 0.05)
#' sims <- redist_flip(map = iowa_map, nsims = 100)
#'
redist_flip <- function(map, nsims, warmup = 0, init_plan,
                        constraints = redist_constr(map) %>% add_constr_edges_rem(0.4),
                        nthin = 1, eprob = 0.05, lambda = 0, temper = FALSE,
                        betaseq = "powerlaw", betaseqlength = 10, betaweights = NULL,
                        adapt_lambda = FALSE, adapt_eprob = FALSE, exact_mh = FALSE,
                        adjswaps = TRUE, init_name = NULL, verbose = TRUE) {
    if (verbose) {
        ## Initialize ##
        cli::cli({
            cli::cli_h1(cli::col_red("redist_flip()"))
            cli::cli_h2(cli::col_red("Automated Redistricting Simulation Using Markov Chain Monte Carlo"))
        })
    }
    # process raw inputs
    nprec <- nrow(map)
    map <- validate_redist_map(map)
    adj <- get_adj(map)
    total_pop <- map[[attr(map, "pop_col")]]
    ndists <- attr(map, "ndists")

    if (any(total_pop >= get_target(map)))
        cli::cli_abort("Units ", which(total_pop >= get_target(map)),
            " have population larger than the district target.\n",
            "Redistricting impossible.")

    # process constraints
    if (!inherits(constraints, "redist_constr")) {
        cli::cli_abort("Not a {.cls redist_constr} object.")
    }


    if (!any(class(nthin) %in% c("numeric", "integer"))) {
        cli::cli_abort("nthin must be an integer")
    } else if (nthin < 1) {
        cli::cli_abort("nthin must be a nonnegative integer.")
    } else {
        nthin <- as.integer(nthin)
    }

    pop_tol <- get_pop_tol(map)


    exist_name <- attr(map, "existing_col")
    if (missing(init_plan)) {
        init_plan <- get_existing(map)

        if (is.null(init_plan)) {
            invisible(capture.output(init_plan <- redist_smc(map,
                nsims = 1,
                silent = TRUE
            ), type = "message"))
            init_plan <- as.matrix(init_plan) - 1L

            if (is.null(init_name)) {
                init_name <- "<init>"
            }

        } else {
            if (is.null(init_name)) init_name <- exist_name
            init_plan <- redist.sink.plan(plan = init_plan)
            components <- contiguity(adj, init_plan)
            if (any(components > 1)) {
                cli::cli_abort("init_plan does not point to a contiguous plan.")
            }
        }
    }

    adapt_lambda <- as.integer(adapt_lambda)
    adapt_eprob <- as.integer(adapt_eprob)
    exact_mh <- as.integer(exact_mh)

    if (verbose) {
        cli::cli_alert_info("Preprocessing data.")
    }

    preprocout <- redist.preproc(
        adj = adj,
        total_pop = total_pop,
        init_plan = init_plan,
        ndists = ndists,
        pop_tol = pop_tol,
        temper = temper,
        betaseq = betaseq,
        betaseqlength = betaseqlength,
        betaweights = betaweights,
        adjswaps = adjswaps,
        maxiterrsg = 1,
        verbose = verbose
    )

    if (verbose) {
        cli::cli_alert_info("Starting swMH().")
    }

    algout <- swMH(
        aList = preprocout$data$adjlist,
        cdvec = preprocout$data$init_plan,
        popvec = preprocout$data$total_pop,
        constraints = as.list(constraints),
        nsims = nsims*nthin + warmup,
        eprob = eprob,
        pct_dist_parity = preprocout$params$pctdistparity,
        beta_sequence = preprocout$params$betaseq,
        beta_weights = preprocout$params$betaweights,
        lambda = lambda,
        beta = preprocout$params$beta,
        adapt_beta = preprocout$params$temperbeta,
        adjswap = preprocout$params$adjswaps,
        exact_mh = exact_mh,
        adapt_lambda = adapt_lambda,
        adapt_eprob = adapt_eprob,
        verbose = as.logical(verbose)
    )

    algout <- redist.warmup.chain(algout, warmup = warmup)
    algout <- redist.thin.chain(algout, thin = nthin)


    algout$plans <- algout$plans + 1

    out <- new_redist_plans(
        plans = algout$plans,
        map = map,
        algorithm = "flip",
        wgt = NULL,
        resampled = FALSE,
        ndists = ndists,
        lambda = lambda,
        eprob = eprob,
        pop_tol = pop_tol,
        adapt_eprob = as.logical(adapt_eprob),
        adapt_lambda = as.logical(adapt_lambda),
        warmup = warmup,
        nthin = nthin,
        mh_acceptance = mean(algout$mhdecisions),
        final_eprob = algout$final_eprob,
        final_lambda = algout$final_lambda,
        distance_original = algout$distance_original,
        beta_sequence = algout$beta_sequence
    ) %>%
        mutate(
            distance_parity = rep(algout$distance_parity, each = ndists),
            mhdecisions = rep(algout$mhdecisions, each = ndists),
            mhprob = rep(algout$mhprob, each = ndists),
            pparam = rep(algout$pparam, each = ndists),
            beta_sequence = rep(algout$beta_sequence, each = ndists),
            energy_psi = rep(algout$energy_psi, each = ndists),
            boundary_partitions = rep(algout$boundary_partitions, each = ndists),
            boundary_ratio = rep(algout$boundary_partitions, each = ndists)
        )
    add_tb <- apply(algout$psi_store, 1, function(x) rep(x, each = ndists)) %>%
        dplyr::as_tibble() %>%
        dplyr::rename_with(function(x) paste0("constraint_", x))

    names_tb <- names(add_tb)[apply(add_tb, 2, function(x) { !all(x == 0) })]
    out <- bind_cols(out, select(add_tb, all_of(names_tb)))

    if (!is.null(init_name) && !isFALSE(init_name)) {
        out <- add_reference(out, init_plan, init_name)
    }

    out
}


#' Flip MCMC Redistricting Simulator using Simulated Annealing
#'
#' \code{redist_flip_anneal} simulates congressional redistricting plans
#' using Markov chain Monte Carlo methods coupled with simulated annealing.
#'
#' @param map A \code{\link{redist_map}} object.
#' @param nsims The number of samples to draw, not including warmup.
#' @param warmup The number of warmup samples to discard.
#' @param init_plan A vector containing the congressional district labels
#' of each geographic unit. The default is \code{NULL}. If not provided,
#' a random initial plan will be generated using \code{redist_smc}. You can also
#' request to initialize using \code{redist.rsg} by supplying 'rsg', though this is
#' not recommended behavior.
#' @param constraints A `redist_constr` object.
#' @param num_hot_steps The number of steps to run the simulator at beta = 0.
#' Default is 40000.
#' @param num_annealing_steps The number of steps to run the simulator with
#' linearly changing beta schedule. Default is 60000
#' @param num_cold_steps The number of steps to run the simulator at beta = 1.
#' Default is 20000.
#' @param eprob The probability of keeping an edge connected. The
#' default is \code{0.05}.
#' @param lambda The parameter determining the number of swaps to attempt
#' each iteration of the algorithm. The number of swaps each iteration is
#' equal to Pois(\code{lambda}) + 1. The default is \code{0}.
#' @param adapt_lambda Whether to adaptively tune the lambda parameter so that the Metropolis-Hastings
#' acceptance probability falls between 20\% and 40\%. Default is FALSE.
#' @param adapt_eprob Whether to adaptively tune the edgecut probability parameter so that the
#' Metropolis-Hastings acceptance probability falls between 20\% and 40\%. Default is
#' FALSE.
#' @param exact_mh Whether to use the approximate (0) or exact (1)
#' Metropolis-Hastings ratio calculation for accept-reject rule. Default is FALSE.
#' @param maxiterrsg Maximum number of iterations for random seed-and-grow
#' algorithm to generate starting values. Default is 5000.
#' @param verbose Whether to print initialization statement.
#' Default is \code{TRUE}.
#'
#' @return redist_plans
#'
#' @concept simulate
#' @export
redist_flip_anneal <- function(map,
                               nsims,
                               warmup = 0,
                               init_plan = NULL,
                               constraints = redist_constr(),
                               num_hot_steps = 40000,
                               num_annealing_steps = 60000,
                               num_cold_steps = 20000,
                               eprob = 0.05,
                               lambda = 0,
                               adapt_lambda = FALSE,
                               adapt_eprob = FALSE,
                               exact_mh = FALSE,
                               maxiterrsg = 5000,
                               verbose = TRUE) {

    if (verbose) {
        ## Initialize ##
        cli::cli({
            cli::cli_h1(cli::col_red("redist_flip_anneal()"))
            cli::cli_h2(cli::col_red("Automated Redistricting Simulation Using Markov Chain Monte Carlo"))
        })
    }
    # process raw inputs
    nprec <- nrow(map)
    map <- validate_redist_map(map)
    adj <- get_adj(map)
    total_pop <- map[[attr(map, "pop_col")]]
    ndists <- attr(map, "ndists")

    if (any(total_pop >= get_target(map))) {
        cli::cli_abort("Units ", which(total_pop >= get_target(map)),
            " have population larger than the district target.\n",
            "Redistricting impossible.")
    }

    # process constraints
    if (!inherits(constraints, "redist_constr")) {
        cli::cli_abort("Not a {.cls redist_constr} object.")
    }


    if (!any(class(nthin) %in% c("numeric", "integer"))) {
        cli::cli_abort("nthin must be an integer")
    } else if (nthin < 1) {
        cli::cli_abort("nthin must be a nonnegative integer.")
    } else {
        nthin <- as.integer(nthin)
    }


    pop_tol <- get_pop_tol(map)


    exist_name <- attr(map, "existing_col")
    if (missing(init_plan)) {
        init_plan <- get_existing(map)

        if (is.null(init_plan)) {
            invisible(capture.output(init_plan <- redist_smc(map,
                nsims = 1,
                silent = TRUE
            ), type = "message"))
            init_plan <- as.matrix(init_plan) - 1L

            if (is.null(init_name)) {
                init_name <- "<init>"
            }

        } else {
            if (is.null(init_name)) init_name <- exist_name
            init_plan <- redist.sink.plan(plan = init_plan)
            components <- contiguity(adj, init_plan)
            if (any(components > 1)) {
                cli::cli_abort("init_plan does not point to a contiguous plan.")
            }
        }
    }

    adapt_lambda <- as.integer(adapt_lambda)
    adapt_eprob <- as.integer(adapt_eprob)
    exact_mh <- as.integer(exact_mh)

    if (verbose) {
        cli::cli_alert_info("Preprocessing data.")
    }

    ## ------------------
    ## Preprocessing data
    ## ------------------
    if (verbose) {
        cat("Preprocessing data.\n\n")
    }
    preprocout <- redist.preproc(adj = adj, total_pop = total_pop,
        init_plan = init_plan, ndists = ndists,
        pop_tol = pop_tol,
        temper = FALSE,
        betaseq = "powerlaw", betaseqlength = 10,
        betaweights = NULL,
        adjswaps = TRUE, maxiterrsg = maxiterrsg,
        verbose = verbose)


    if (verbose) {
        cat("Starting swMH().\n")
    }

    algout <- swMH(aList = preprocout$data$adjlist,
        cdvec = preprocout$data$init_plan,
        popvec = preprocout$data$total_pop,
        constraints = as.list(constraints),
        nsims = 100,
        eprob = eprob,
        pct_dist_parity = preprocout$params$pctdistparity,
        beta_sequence = preprocout$params$betaseq,
        beta_weights = preprocout$params$betaweights,
        lambda = lambda,
        beta = 0,
        adapt_beta = "annealing",
        adjswap = preprocout$params$adjswaps,
        exact_mh = exact_mh,
        adapt_lambda = adapt_lambda,
        adapt_eprob = adapt_eprob,
        num_hot_steps = num_hot_steps,
        num_annealing_steps = num_annealing_steps,
        num_cold_steps = num_cold_steps,
        verbose = as.logical(verbose))

    algout$plans <- algout$plans + 1


    out <- new_redist_plans(
        plans = algout$plans,
        map = map,
        algorithm = "flip_anneal",
        wgt = NULL,
        resampled = FALSE,
        ndists = ndists,
        lambda = lambda,
        eprob = eprob,
        pop_tol = pop_tol,
        adapt_eprob = as.logical(adapt_eprob),
        adapt_lambda = as.logical(adapt_lambda),
        warmup = warmup,
        nthin = nthin,
        mh_acceptance = mean(algout$mhdecisions)
    ) %>% mutate(
        distance_parity = rep(algout$distance_parity, each = ndists),
        mhdecisions = rep(algout$mhdecisions, each = ndists),
        mhprob = rep(algout$mhprob, each = ndists),
        pparam = rep(algout$pparam, each = ndists),
        beta_sequence = rep(algout$beta_sequence, each = ndists),
        energy_psi = rep(algout$energy_psi, each = ndists),
        boundary_partitions = rep(algout$boundary_partitions, each = ndists),
        boundary_ratio = rep(algout$boundary_partitions, each = ndists)
    )
    add_tb <- apply(algout$psi_store, 1, function(x) rep(x, each = ndists)) %>%
        dplyr::as_tibble() %>%
        dplyr::rename_with(function(x) paste0("constraint_", x))

    names_tb <- names(add_tb)[apply(add_tb, 2, function(x) { !all(x == 0) })]
    out <- bind_cols(out, select(add_tb, all_of(names_tb)))

    if (!is.null(init_name) && !isFALSE(init_name)) {
        out <- add_reference(out, init_plan, init_name)
    }

    out
}
