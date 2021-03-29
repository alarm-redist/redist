#' Static Variation of Information Plot
#'
#' @param plans matrix of district assignments
#' @param group_pop Required Population of subgroup being studied in each precinct.
#' @param total_pop Required. Population of each precinct.
#' @param shp sf dataframe
#'
#'
#' @importFrom  patchwork wrap_plots plot_layout
#' @importFrom dplyr summarize
#' @importFrom ggplot2 geom_text theme_bw
#' @importFrom stats cmdscale kmeans
#'
#' @return patchworked ggplot
#' @concept plot
#' @export
redist.plot.varinfo <- function(plans, group_pop, total_pop, shp){
  centers <- 5
  gp <- redist.group.percent(plans = plans, group_pop = group_pop, total_pop = total_pop)
  pct_min <- colmin(gp)
  dists <- redist.distances(plans, 'info', total_pop=total_pop)$VI
  mds <- cmdscale(dists)
  tb <- tibble(mds1 = mds[,1], mds2 = mds[,2],
               pct_min = pct_min, id = 1:ncol(plans))

  km <- kmeans(x = mds, centers = centers)

  tb$cluster <- km$cluster

  tb$dist <- NA_real_
  for(i in 1:nrow(tb)){
    tb$dist[i] <- sqrt( (tb$mds1[i] - km$centers[tb$cluster[i], 1])^2 + (tb$mds2[i] - km$centers[tb$cluster[i], 2])^2  )
  }

  sub <- tb %>% group_by(cluster) %>% filter(dist == min(dist)) %>% slice(1)

  p <- tb %>%
    ggplot(aes(x = mds1, y = mds2, color = pct_min)) +
    geom_point() +
    theme_bw() +
    labs(x = '', y = '', color = 'Minority\nPercent') +
    geom_text(data = sub, aes(label = id), color = 'red') +
    ggplot2::scale_color_viridis_c()


  ratio <- group_pop/total_pop

  p2  <- lapply(1:(centers+1), function(x){
    if(x == 1){
      return(NULL)
    }
    shp$newcd  <- as.character(plans[,sub$id[x-1]])
    shpdist <- shp %>% group_by(newcd) %>% summarize(geometry = st_union(geometry))
    shp %>% ggplot() +
      geom_sf(aes(fill = ratio)) +
      labs(fill = 'Minority\nPercent', title = sub$id[x-1]) +
      theme_void() +
      scale_fill_gradient(limits = c(0,1), low = '#ffffff', high = '#08306b') +
      geom_sf(data = shpdist, fill = NA, lwd = 1, color = 'orange')
  })

  p2[[1]] <- p

  layouts <- "
  AAB
  AAC
  DEF
  "

  return(patchwork::wrap_plots(p2) + patchwork::plot_layout(design = layouts, guides = 'collect'))
}

utils::globalVariables(c('geometry', 'cluster', 'mds1', 'mds2', 'id', 'newcd'))





#' Make a comparison between two sets of plans
#'
#' This function provides one way to identify the structural differences between
#' two sets of redistricting plans. It operates by computing the precinct
#' co-occurrence matrix (a symmetric matrix where the i,j-th entry is the
#' fraction of plans where precinct i and j are in the same district) for each
#' set, and then computing the first eigenvalue of the difference in these two
#' matrices (in each direction). These eigenvalues identify the important parts
#' of the map.
#'
#' The co-occurence matrices are regularized with a \eqn{Beta(1/ndists, 1-1/ndists)}
#' prior, which is useful for when either `set1` or `set2` is small.
#'
#' @param plans a [redist_plans] object
#' @param set1,set2 [`<data-masking>`][dplyr::dplyr_data_masking] indexing vectors
#'   for the plan draws to compare. Must be mutually exclusive.
#' @param thresh the value to threshold the eigenvector at in determining the
#'   relevant set of precincts for comparison.
#'
#' @return A list with the following elements:
#' \item{eigen1}{A numeric vector containing the first eigenvector of
#'   \code{p1 - p2}, where \code{p1} and \code{p2} are the co-occurrence matrices
#'   for \code{set1} and \code{set2}, respectively.}
#' \item{eigen2}{A numeric vector containing the first eigenvector of
#'   \code{p2 - p1}, where \code{p1} and \code{p2} are the co-occurrence matrices
#'   for \code{set1} and \code{set2}, respectively.}
#' \item{group_1a, group_1b}{Lists of precincts. Compared to \code{set2}, in the
#'   \code{set1} plans these precincts were much more likely to be in separate
#'   districts. Computed by thresholding \code{eigen1} at \code{thresh}.}
#' \item{group_2a, group_2b}{Lists of precincts. Compared to \code{set1}, in the
#'   \code{set2} plans these precincts were much more likely to be in separate
#'   districts. Computed by thresholding \code{eigen2} at \code{thresh}.}
#' \item{cooccur_sep_1}{The difference in the average co-occurrence of precincts
#'   in \code{group_1a} and \code{group_1b} between \code{set2} and \code{set1}.
#'   Higher indicates better separation.}
#' \item{cooccur_sep_2}{The difference in the average co-occurrence of precincts
#'   in \code{group_2a} and \code{group_2b} between \code{set1} and \code{set2}.
#'   Higher indicates better separation.}
#'
#'
#' @md
#' @concept analyze
#' @export
compare_plans = function(plans, set1, set2, thresh=0.1) {
    stopifnot(inherits(plans, "redist_plans"))

    set1 = eval_tidy(enquo(set1), plans)
    set2 = eval_tidy(enquo(set2), plans)
    if (is.logical(set1)) set1 = which(set1)
    if (is.logical(set2)) set2 = which(set2)
    if (length(intersect(set1, set2)) > 0)
        stop("`set1` and `set2` must be mutually exclusive.")
    n1 = length(set1)
    n2 = length(set2)

    pm = get_plan_matrix(plans)
    base_co = 1 / max(pm[, 1]) # baseline coccurence

    p1 = (n1*prec_cooccur(pm, set1) + base_co) / (n1 + 1)
    p2 = (n2*prec_cooccur(pm, set2) + base_co) / (n2 + 1)

    evec1 = eigen(p1 - p2, symmetric=TRUE)$vectors[, 1]
    evec2 = eigen(p2 - p1, symmetric=TRUE)$vectors[, 1]

    group_1a = which(evec1 >= thresh)
    group_1b = which(evec1 <= -thresh)
    group_2a = which(evec2 >= thresh)
    group_2b = which(evec2 <= -thresh)

    list(eigen1 = evec1,
         eigen2 = evec2,
         group_1a = group_1a,
         group_1b = group_1b,
         group_2a = group_2a,
         group_2b = group_2b,
         cooccur_sep_1 = mean(p2[group_1a, group_1b]) -
             mean(p1[group_1a, group_1b]),
         cooccur_sep_2 = mean(p1[group_2a, group_2b]) -
             mean(p2[group_2a, group_2b]))
}

