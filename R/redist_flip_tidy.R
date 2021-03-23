#' 'Flip' Markov Chain Monte Carlo Redistricting Simulation
#'
#' \code{redist_flip} provides an improved interface to the methods in
#' \code{\link{redist.mcmc}}.
#'
#' This function allows users to simulate redistricting plans
#' using the Markov Chain Monte Carlo methods of Fifield et al. Several
#' constraints corresponding to substantive requirements in the redistricting
#' process are implemented, including population parity and geographic
#' compactness. In addition, the function includes multiple-swap and simulated
#' tempering functionality to improve the mixing of the Markov Chain.
#'
#' @param map A \code{\link{redist_map}} object.
#' @param nsims The number of samples to draw.
#' @param init_plan A vector containing the congressional district labels
#' of each geographic unit. The default is \code{NULL}. If not provided,
#' a random initial plan will be generated using \code{redist.smc}.
#' @param counties A column in map containing county membership
#' @param group_pop A column in map containing group populations
#' @param constraints a list of constraints to implement. Can be created with
#' \code{flip_constraints_helper}
#' @param eprob The probability of keeping an edge connected. The
#' default is \code{0.05}.
#' @param lambda lambda The parameter determining the number of swaps to attempt
#' each iteration of the algorihtm. The number of swaps each iteration is
#' equal to Pois(\code{lambda}) + 1. The default is \code{0}.
#' @param temper Whether to use simulated tempering algorithm. Default is FALSE.
#' @param betaseq Sequence of beta values for tempering. The default is
#' \code{powerlaw} (see Fifield et. al (2020) for details).
#' @param betaseqlength @param betaseqlength Length of beta sequence desired for
#' tempering. The default is \code{10}.
#' @param betaweights betaweights Sequence of weights for different values of
#' beta. Allows the user to upweight certain values of beta over
#' others. The default is \code{NULL} (equal weighting).
#' @param adapt_lambda adapt_lambda Whether to adaptively tune the lambda parameter so that the Metropolis-Hastings
#' acceptance probability falls between 20\% and 40\%. Default is FALSE.
#' @param adapt_eprob eprob Whether to adaptively tune the edgecut probability parameter so that the
#' Metropolis-Hastings acceptance probability falls between 20\% and 40\%. Default is
#' FALSE.
#' @param exact_mh  Whether to use the approximate (FALSE) or exact (TRUE)
#' Metropolis-Hastings ratio calculation for accept-reject rule. Default is FALSE.
#' @param adjswaps Flag to restrict swaps of beta so that only
#' values adjacent to current constraint are proposed. The default is
#' \code{TRUE}.
#' @param init_name a name for the initial plan, or \code{FALSE} to not include
#'   the initial plan in the output.  Defaults to the column name of the
#'   existing plan, or "\code{<init>}" if the initial plan is sampled.
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
#' @importFrom rlang eval_tidy enquo
#' @importFrom utils capture.output
#'
#' @examples \dontrun{
#' data(iowa)
#' iowa_map <- redist_map(iowa, ndists = 4, existing_plan = cd_2010, total_pop = 'pop')
#' sims <- redist_flip(map = iowa_map, nsims = 100)
#' }
redist_flip <- function(map, nsims, init_plan, counties = NULL, group_pop, constraints = list(),
                        eprob = 0.05, lambda = 0, temper = FALSE,
                        betaseq = 'powerlaw', betaseqlength = 10, betaweights = NULL,
                        adapt_lambda = FALSE, adapt_eprob = FALSE, exact_mh = FALSE,
                        adjswaps = TRUE, init_name = NULL, verbose = TRUE) {
  if (verbose) {
    ## Initialize ##
    divider <- c(paste(rep('=', 20), sep = '', collapse = ''), '\n')

    cat('\n', append = TRUE)
    cat(divider, append = TRUE)
    cat('redist.mcmc(): Automated Redistricting Simulation Using
         Markov Chain Monte Carlo\n\n', append = TRUE)
  }
  nprec <- nrow(map)
  map <- validate_redist_map(map)
  adj <- get_adj(map)
  total_pop <- map[[attr(map, 'pop_col')]]
  ndists <- attr(map, 'ndists')

  if ('pop_tol' %in% names(attributes(map))) {
    pop_tol <- attr(map, 'pop_tol')
  } else if ('pop_bounds' %in% names(attributes(map))) {
    pop_tol <- attr(map, 'pop_bounds')[3] / attr(map, 'pop_bounds')[2] - 1
  } else {
    pop_tol <- 100
  }

  counties <- eval_tidy(enquo(counties), map)
  if (is.null(counties)) {
    counties <- rep(1, length(adj))
  } else {
    counties <- redist.county.id(counties)
    components <- contiguity(adj, counties)
    if (any(components > 1)) {
      warning('counties are not contiguous. Expect additional splits.')
      counties <- redist.county.id(redist.county.relabel(adj, counties))
    }
  }

  exist_name = attr(map, "existing_col")
  if (missing(init_plan)) {
    init_plan <- get_existing(map)

    if (is.null(init_plan)) {
      invisible(capture.output(init_plan <- redist.smc(
        adj = adj,
        total_pop = total_pop,
        nsims = 1,
        ndists = ndists,
        counties = counties,
        pop_tol = pop_tol,
        silent = TRUE
      ), type = 'message'))
      init_plan <- init_plan$plans
      if (is.null(init_name)) init_name = "<init>"
    } else {
      if (is.null(init_name)) init_name = exist_name
      init_plan <- redist.sink.plan(plan = init_plan)
      components <- contiguity(adj, init_plan)
      if (any(components > 1)) {
        stop('init_plan does not point to a contiguous plan.')
      }
    }
  }

  if (missing(group_pop)) {
    group_pop <- rep(0, nprec)
  } else {
    group_pop <- eval_tidy(enquo(group_pop), map)
  }

  pre_pre_proc <- process_flip_constr(constraints, group_pop, counties)

  constraints_name <- names(pre_pre_proc)
  constraints_wt <- sapply(pre_pre_proc, function(x) {
    x$weight
  })
  if (all(constraints_wt == 0)) {
    constraints_name <- constraints_wt <- NULL
  } else {
    constraints_name <- constraints_name[constraints_wt != 0]
    constraints_wt <- constraints_wt[constraints_wt != 0]
  }


  adapt_lambda <- as.integer(adapt_lambda)
  adapt_eprob <- as.integer(adapt_eprob)
  exact_mh <- as.integer(exact_mh)

  if (verbose) {
    cat('Preprocessing data.\n\n')
  }
  preprocout <- redist.preproc(
    adj = adj,
    total_pop = total_pop,
    init_plan = init_plan,
    ndists = ndists,
    pop_tol = pop_tol,
    counties = counties,
    group_pop = group_pop,
    areasvec = pre_pre_proc$compact$areasvec,
    borderlength_mat = pre_pre_proc$compact$borderlength_mat,
    ssdmat = pre_pre_proc$compact$ssdmat,
    compactness_metric = pre_pre_proc$compact$metric,
    partisan_metric = pre_pre_proc$partisan$metric,
    temper = temper,
    constraint = constraints_name,
    constraintweights = constraints_wt,
    betaseq = betaseq,
    betaseqlength = betaseqlength,
    betaweights = betaweights,
    adjswaps = adjswaps,
    maxiterrsg = 1,
    contiguitymap = 'rooks',
    tgt_min = pre_pre_proc$vra$target_min,
    tgt_other = pre_pre_proc$vra$target_other,
    rvote = pre_pre_proc$partisan$rvote,
    dvote = pre_pre_proc$partisan$dvote,
    minorityprop = pre_pre_proc$hinge$minorityprop
  )

  if (verbose) {
    cat('Starting swMH().\n')
  }
  algout <- swMH(
    aList = preprocout$data$adjlist,
    cdvec = preprocout$data$init_plan,
    cdorigvec = preprocout$data$init_plan,
    popvec = preprocout$data$total_pop,
    grouppopvec = preprocout$data$group_pop,
    areas_vec = preprocout$data$areasvec,
    county_membership = preprocout$data$counties,
    borderlength_mat = preprocout$data$borderlength_mat,
    nsims = nsims,
    eprob = eprob,
    pct_dist_parity = preprocout$params$pctdistparity,
    beta_sequence = preprocout$params$betaseq,
    beta_weights = preprocout$params$betaweights,
    ssdmat = preprocout$data$ssdmat,
    lambda = lambda,
    beta = preprocout$params$beta,
    weight_population = preprocout$params$weightpop,
    weight_compact = preprocout$params$weightcompact,
    weight_segregation = preprocout$params$weightseg,
    weight_vra = preprocout$params$weightvra,
    weight_similar = preprocout$params$weightsimilar,
    weight_countysplit = preprocout$params$weightcountysplit,
    weight_partisan = preprocout$params$weightpartisan,
    weight_minority = preprocout$params$weightminority,
    weight_hinge = preprocout$params$weighthinge,
    adapt_beta = preprocout$params$temperbeta,
    adjswap = preprocout$params$adjswaps,
    exact_mh = exact_mh,
    adapt_lambda = adapt_lambda,
    adapt_eprob = adapt_eprob,
    compactness_measure = pre_pre_proc$compact$metric,
    partisan_measure = preprocout$params$partisan_metric,
    ssd_denom = pre_pre_proc$compact$ssd_denom,
    tgt_min = pre_pre_proc$vra$target_min,
    tgt_other = pre_pre_proc$vra$target_other,
    rvote = preprocout$params$rvote,
    dvote = preprocout$params$dvote,
    minorityprop = preprocout$params$minorityprop,
    verbose = as.logical(verbose)
  )


  if (min(algout$plans) == 0) {
    algout$plans <- ifelse(algout$plans == 0, ndists, algout$plans)
  }

  out <- new_redist_plans(
    plans = algout$plans,
    map = map,
    algorithm = 'mcmc',
    wgt = NULL,
    resampled = NULL,
    lambda = lambda,
    eprob = eprob,
    adapt_eprob = as.logical(adapt_eprob),
    adapt_lambda = as.logical(adapt_lambda)
  ) %>% mutate(
    boundary_partitions = rep(algout$boundary_partitions, each = ndists),
    boundary_ratio = rep(algout$boundary_partitions, each = ndists)
  )
  add_tb <- tibble(
    constraint_population = rep(algout$constraint_pop, each = ndists),
    constraint_compact = rep(algout$constraint_compact, each = ndists),
    constraint_segregation = rep(algout$constraint_segregation, each = ndists),
    constraint_vra = rep(algout$constraint_vra, each = ndists),
    constraint_similar = rep(algout$constraint_similar, each = ndists),
    constraint_countysplit = rep(algout$constraint_countysplit, each = ndists),
    constraint_partisan = rep(algout$constraint_partisan, each = ndists),
    constraint_minority = rep(algout$constraint_minority, each = ndists),
    constraint_hinge = rep(algout$constraint_hinge, each = ndists)
  )

  names_tb = names(add_tb)[apply(add_tb, 2, function(x) !all(x == 0) )]
  out <- dplyr::bind_cols(out, select(add_tb, all_of(names_tb)))

  if (!is.null(init_name) && !isFALSE(init_name)) {
      out <- add_reference(out, init_plan, init_name)
    }

  return(out)
}


# Helper constraint processor -- flip
# Constraint defaults contained elsewhere but also kinda here
#
# @param constraints passed into `redist_flip`
#
# @return a list with new `constraints`
process_flip_constr <- function(constraints, group_pop, counties) {
  defaults <- list(
    compact = list(
      weight = 0.6, metric = 'edges-removed',
      areas = c(0, 0, 0, 0), borderlength_mat = matrix(0, 2, 2),
      ssdmat = matrix(1, 2, 2), ssd_denom = 1.0
    ),
    population = list(weight = 0),
    countysplit = list(weight = 0, counties = counties),
    hinge = list(weight = 0, minorityprop = 0.55, group_pop = group_pop),
    vra = list(
      weight = 0, target_min = 0.55, target_other = 0.25,
      pow_vra = 1.5, group_pop = group_pop
    ),
    minority = list(weight = 0, minorityprop = 0.55, group_pop = group_pop),
    similarity = list(weight = 0),
    partisan = list(weight = 0, rvote = 0, dvote = 0, metric = 'efficiency-gap'),
    segregation = list(weight = 0, group_pop = group_pop)
  )

  for (type in names(constraints)) {
    for (el in names(constraints[[type]])) {
      defaults[[type]][[el]] <- constraints[[type]][[el]]
    }
  }

  # if(defaults$compact$weight > 0 & defaults$compact$metric == 'polsby-popper'){
  #   if(all(defaults$compact$areas == 0)) {
  #     stop('constraint compact used with polsby-popper without areas entry.')
  #   }
  #   if(all(defaults$compact$areas == 0)) {
  #     stop('constraint compact used with polsby-popper without areas entry.')
  #   }
  # }

  return(defaults)
}
