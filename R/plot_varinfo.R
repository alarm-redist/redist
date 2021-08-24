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
    shpdist <- shp %>% group_by(newcd) %>% summarize(geometry = st_union(sf::st_geometry(geometry)))
    shp %>% ggplot() +
      geom_sf(aes(fill = ratio)) +
      labs(fill = 'Minority\nPercent', title = sub$id[x-1]) +
      theme_void() +
      ggplot2::scale_fill_gradient(limits = c(0,1), low = '#ffffff', high = '#08306b') +
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



