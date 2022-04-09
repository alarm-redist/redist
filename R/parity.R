#' Calculates Maximum Deviation from Population Parity
#'
#' Computes the deviation from population parity from a plan.
#' Higher values indicate that (at least) a single district in the map deviates
#' from population parity. See Details.
#'
#' @details With a map with \code{pop} representing the populations of each district,
#'  the deviation from population parity is given as \code{max(abs(pop - parity) / parity)}
#'  where \code{parity = sum(pop)/length(pop)} is the population size for the
#'  average district.
#'  Therefore, the metric can be thought of as the maximum percent deviation from
#'  equal population. For example, a value of 0.03 in this metric indicates that
#'  all districts are within 3 percent of population parity.
#'
#' @param plans A matrix with one row for each precinct and one column for each
#'   map. Required.
#' @param total_pop A numeric vector with the population for every precinct.
#' @param ncores Number of cores to use for parallel computing. Default is 1.
#'
#' @importFrom foreach %do% %dopar% foreach
#' @return numeric vector with the population parity for each column
#'
#' @concept analyze
#' @export
redist.parity <- function(plans, total_pop, ncores = 1) {
    if (!is.numeric(total_pop)) {
        cli_abort("{.arg total_pop} must be a numeric vector")
    }
    if (!is.matrix(plans)) {
        plans <- matrix(plans, ncol=1)
    }
    if (!is.matrix(plans)) {
        cli_abort("{.arg plans} must be a matrix")
    }

    if (length(total_pop) != nrow(plans)) {
        cli_abort(".arg plans} and {.arg total_pop} must have same number of precincts.")
    }

    # parallelize as in fastLink package to avoid Windows/unix issues
    N = ncol(plans)
    nc <- min(ncores, max(1, floor(N/2)))

    if (nc == 1){
        `%oper%` <- `%do%`
    } else {
        `%oper%` <- `%dopar%`
        cl <- makeCluster(nc, setup_strategy = 'sequential', methods=FALSE)
        registerDoParallel(cl)
        on.exit(stopCluster(cl))
    }

    if (min(plans[,1]) == 0)
        plans = plans + 1
    n_distr = max(plans[,1])


    chunks = split(1:N, rep(1:nc, each=ceiling(N/nc))[1:N])
    out = foreach(map=chunks, .combine = "c") %oper% {
        max_dev(plans[, map, drop = FALSE], total_pop, n_distr)
    }

    unlist(out)
}


#' Calculates Sparse Population Moves to Minimize Population Deviation
#'
#' This function computes a minimal set of population moves (e.g., 5 people from
#' district 1 to district 3) to maximally balance the population between
#' districts. The moves are only allowed between districts that share the
#' territory of a county, so that any boundary adjustments are guaranteed to
#' preserve all unbroken county boundaries.
#'
#' @param map a [redist_map]
#' @param plan an integer vector containing the plan to be balanced.
#'   Tidy-evaluated.
#' @param counties an optional vector of counties, whose boundaries will be
#'   preserved. Tidy-evaluated.
#' @param penalty the larger this value, the more to encourage sparsity.
#'
#' @returns a list with components:
#' \describe{
#'   \item{`moves`}{A tibble describing the population moves}
#'   \item{`pop_old`}{The current district populations}
#'   \item{`pop_new`}{The district populations after the moves}
#' }
#'
#' @concept analyze
#' @md
#' @export
min_move_parity = function(map, plan, counties=NULL, penalty=0.2) {
    adj = get_adj(map)
    V = length(adj)
    nd = attr(map, "ndists")

    plan = eval_tidy(enquo(plan), map)
    if (!is.numeric(plan) && all(plan > 0) && length(plan) == V)
        cli_abort("{.arg plan} must be a positive integer vector with one entry per precinct.")

    if (missing(counties)) {
        counties = rep(1L, V)
    } else {
        counties = as.integer(as.factor(eval_tidy(enquo(counties), map)))
    }

    distr_adj = get_plan_graph(adj, length(adj), plan, nd)
    edges = do.call(rbind, lapply(seq_along(distr_adj), function(i) {
        tibble(from=i, to=distr_adj[[i]] + 1L)
    })) %>%
        rowwise() %>%
        filter(.data$from < .data$to,
               any(unique(counties[plan == .data$from]) %in% counties[plan == .data$to])) %>%
        ungroup()

    n_edge = nrow(edges)
    e_idx = as.matrix(mutate(edges, i = row_number()))
    diff_mat = matrix(0, nrow=nd, ncol=n_edge)
    diff_mat[e_idx[, -2]] = -1
    diff_mat[e_idx[, -1]] = 1

    pops = pop_tally(matrix(plan, ncol=1), map[[attr(map, "pop_col")]], nd)
    discrep = pops - mean(pops)

    fn_balance = function(move, alpha=0.1) {
        sum(abs(discrep + diff_mat %*% move)) + alpha*sum(abs(move))
    }
    gr_balance = function(move, alpha=0.1) {
        t(sign(discrep + diff_mat %*% move)) %*% diff_mat + alpha * sign(move)
    }

    res = optim(rep(0, n_edge), fn_balance, gr_balance, alpha=penalty,
                method="BFGS", control=list(maxit=1e3, reltol=1e-9, abstol=1))

    move = round(res$par)
    pops_new = pops + diff_mat %*% move
    from_old = edges$from
    edges = mutate(edges,
                   from = if_else(move < 0, .data$to, .data$from),
                   to = if_else(move < 0, from_old, .data$to),
                   move = abs(move)) %>%
        filter(.data$move > 0)

    list(moves = edges,
         pop_old = pops[,1],
         pop_new = pops_new[,1])
}
