#' Pick One Plan from Many Plans
#'
#' @param plans a redist_plans object
#' @param map a redist_map object
#' @param counties A column in map with county names. Defaults to NULL and will
#' assume each row in map is its own county if left NULL
#' @param comp A column in plans with compactness. Defaults to NULL and will
#' calculate FracKept if left NULL
#' @param maximize_comp Should comp be maximized? Defaults to TRUE.
#' @param cut_point portion of plans to keep at each step
#'
#' @return An integer vector with a single plan
#' @export
#'
#' @importFrom stats quantile
#'
#' @examples
#' data(iowa)
#' map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01, total_pop = pop)
#' sims <- redist_smc(map, 100)
#' pick_a_plan(sims, map)
#'
pick_a_plan <- function(plans, map, counties = NULL, comp = NULL,
                        maximize_comp = TRUE,  cut_point = 0.3){
  if (!is.numeric(cut_point)) {
    stop('cut_point must be numeric')
  }
  if (cut_point <=0 | cut_point >= 1) {
    stop('cut_point must be greater than zero and less than 1.')
  }

  plans <- validate_redist_plans(x = plans)
  plans <- subset_sampled(plans = plans)
  map <- validate_redist_map(data = map)

  comp <- eval_tidy(enquo(comp), map)
  if (is.null(comp)) {
    comp <- redist.compactness(adj = get_adj(map), plans = get_plans_matrix(plans),
                               measure = 'FracKept')$FracKept
  }

  plans <- plans %>% mutate(compactness = comp,
                            parity = plan_parity(map))

  counties <- eval_tidy(enquo(counties), map)
  if (!is.null(counties)) {
    plans <- plans %>% mutate(splits = county_splits(map, counties))
  } else {
    plans <- plans %>% mutate(splits = 0)
  }

  orig_rows <- prior_rows <- nrow(plans)
  while (prior_rows > 0.05 * orig_rows) {

    plans <- plans %>% filter(.data$parity <= quantile(.data$parity, cut_point))

    if (maximize_comp) {
      plans <- plans %>% filter(.data$compactness >= quantile(comp, 1 - cut_point))
    } else {
      plans <- plans %>% filter(.data$compactness <= quantile(comp, cut_point))
    }

    plans <- plans %>% filter(splits <= quantile(splits, cut_point))


    if (prior_rows == nrow(plans)) break
    prior_rows <- nrow(plans)
  }

  best_draw <- as.integer(plans$draw[which.min(plans$parity)])

  return(get_plans_matrix(plans)[, best_draw])
}


#' Local Plan Optimization
#'
#' Searches the local area for a combination of minimizing county splits, compactness,
#' population parity, and keeping close to the original plan
#'
#' @param plan a single plan to optimize from
#' @param map a redist map object
#' @param counties Required
#'
#' @return a redist_plans object with one plan
#' @export
#'
#' @examples
#' \donttest{
#' data(iowa)
#' map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01, total_pop = pop)
#' plan <- get_plans_matrix(redist_smc(map, 1))[,2]
#' local <- persily(plan = plan, map = map, counties = region)
#' }
persily <- function(plan, map, counties = NULL){

  counties <- eval_tidy(enquo(counties), map)
  cons <- flip_constraints_helper(map = map, init_plan = plan,
                                  constraint = c('compact', 'similarity', 'population', 'countysplit'),
                                  counties = counties,
                                  constraintweight = c(0.1, 10, 10, 10))

  par <- redist.parity(plan, total_pop = map[[attr(map, 'pop_col')]])
  map <- set_pop_tol(map, 2*par)

  bursts <- redist_shortburst(map, score_fn = (scorer_frac_kept(map = map) - scorer_pop_dev(map = map) -
                                                 scorer_splits(map = map, counties = counties) +
                                                 scorer_status_quo(map = map, existing_plan = plan)),
                              backend = 'flip',
                              max_bursts = 100, flip_constraints = cons, return_all = FALSE)
  return(bursts)
}


#coalesce_plans <- function(plan1, plan2, ref_plan = FALSE){
#  return('sorry not yet implemented')
#}
