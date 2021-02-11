#' Static Variation of Information Plot
#'
#' @param district_membership matrix of district assignments
#' @param grouppop Required. Population of subgroup being studied in each precinct.
#' @param fullpop Required. Population of each precinct.
#' @param shp sf dataframe
#' @importFrom  patchwork wrap_plots plot_layout
#' @importFrom viridis scale_color_viridis
#' @importFrom dplyr summarize
#' @importFrom ggplot2 geom_text theme_bw
#' @importFrom stats cmdscale kmeans
#' 
#' @return patchworked ggplot
#' @export
redist.varinfo.plot <- function(district_membership, grouppop, fullpop, shp){
  centers <- 5
  gp <- redist.group.percent(district_membership = district_membership, grouppop = grouppop, fullpop = fullpop)
  pct_min <- apply(gp, 2, function(x){min(x)})
  dists <- redist.distances(district_membership, 'info', pop=grouppop)$VI
  mds <- cmdscale(dists)
  tb <- tibble(mds1 = mds[,1], mds2 = mds[,2], 
               pct_min = pct_min, id = 1:ncol(district_membership)) 
  
  
  
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
    viridis::scale_color_viridis()
  
  
  ratio <- grouppop/fullpop
  
  p2  <- lapply(1:(centers+1), function(x){
    if(x == 1){
      return(NULL)
    }
    shp$newcd  <- as.character(district_membership[,sub$id[x-1]])
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

globalVariables(c('geometry', 'cluster', 'mds1', 'mds2', 'id', 'newcd'))