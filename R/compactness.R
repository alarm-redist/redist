##############################################
## Author: Christopher T Kenny
## Institution: Harvard University
## Date Created: 2020/01/20
## Date Modified: 2020/08/20
## Purpose: R function to compute compactness
##############################################


#' Calculate compactness measures for a set of districts
#' 
#' \code{redist.compactness} is used to compute different compactness statistics for a 
#' shapefile. It currently computes the Polsby-Popper, Schwartzberg score, Length-Width Ratio,
#' Convex Hull score, Reock score, Boyce Clark Index, Fryer Holden score, Edges Removed number, 
#' and the log of the Spanning Trees.
#' 
#' @param shp A SpatialPolygonsDataFrame or sf object. Required unless "EdgesRemoved"
#' and "logSpanningTree" with adjacency provided.
#' @param district_membership A numeric vector (if only one map) or matrix with one row 
#' for each precinct and one column for each map. Required.
#' @param measure A vector with a string for each measure desired. "PolsbyPopper", 
#' "Schwartzberg", "LengthWidth", "ConvexHull", "Reock", "BoyceClark", "FryerHolden",  
#' "EdgesRemoved", and "logSpanningTree" are implemented. Defaults to "PolsbyPopper". Use "all" to
#' return all implemented measures.
#' @param population A numeric vector with the population for every observation. Is
#' only necessary when "FryerHolden" is used for measure. Defaults to NULL.
#' @param adjacency Recommended. A zero-indexed adjacency list. Only used for "PolsbyPopper",
#' EdgesRemoved" and "logSpanningTree". Created with \code{redist.adjacency} if not 
#' supplied and needed. Default is NULL. 
#' @param nloop A numeric to specify loop number. Defaults to 1 if only one map provided 
#' and the column number if multiple maps given.
#' @param ncores Number of cores to use for parallel computing. Default is 1.
#' @param counties A numeric vector from 1:ncounties corresponding to counties. Required for "logSpanningTree".
#' @param pprcpp Boolean, whether to run Polsby Popper using Rcpp. It has a high upfront cost, but quickly becomes faster. 
#' Becomes TRUE if ncol(district_membership > 8).
#' @details This function computes specified compactness scores for a map.  If 
#' there is more than one shape specified for a single district, it combines 
#' them, if necessary, and computes one score for each district.
#' 
#' Polsby-Popper is computed as \deqn{\frac{4*\pi*A(d)}{P(d)^2}} where A is the area 
#' function, the district is d, and P is the perimeter function.
#' 
#' Schwartzberg is computed as \deqn{\frac{P(d)}{2*\pi*\sqrt{\frac{A(d)}{\pi}}}}
#' where A is the area function, the district is d, and P is the perimeter function.
#' 
#' The Length Width ratio is computed as \deqn{\frac{length}{width}} where length
#' is the shorter of the maximum x distance and the maximum y distance. Width is 
#' the longer of the two values.
#' 
#' The Reock score is computed as \deqn{\frac{A(d)}{A(CVH)}} where A is the area
#' function, d is the district, and CVH is the convex hull of the district.
#' 
#' The Boyce Clark Index is computed as \deqn{1 - \sum_{1}^{16}\{\frac{|\frac{r_i}{\sum_ir_i}*100-6.25 |\}}{200}}.
#' The \eqn{r_i} are the distances of the 16 radii computed from the geometric 
#' centroid of the shape to the most outward point of the shape that intersects
#' the radii, if the centroid is contained within the shape.  If the centroid
#' lies outside of the shape, a point on the surface is used, which will naturally
#' incur a penalty to the score.
#' 
#' The Fryer Holden score for each district is computed with \deqn{Pop\odot D(precinct)^2},
#' where \eqn{Pop} is the population product matrix.  Each element is the 
#' product of the ith and jth precinct's populations.  D represents the distance, 
#' where the matrix is the distance between each precinct.  To fully compute this 
#' index, for any map, the sum of these values should be used as the numerator. 
#' The denominator can be calculated from the full enumeration of districts as the
#' smallest calculated numerator.
#' 
#' The log spanning tree measure is the log number of spanning trees.
#' 
#' The edges removed measure is number of egdes removed from the underlying adjacency graph.
#' 
#' @return A tibble with a column that specifies the district, a column for 
#' each specified measure, and a column that specifies the map number.
#' 
#' @references Boyce, R., & Clark, W. 1964. The Concept of Shape in Geography. 
#' Geographical Review, 54(4), 561-572.
#' 
#' Cox, E. 1927. A Method of Assigning Numerical and Percentage Values to the 
#' Degree of Roundness of Sand Grains. Journal of Paleontology, 1(3), 179-183.
#' 
#' Fryer R, Holden R. 2011. Measuring the Compactness of Political Districting Plans. 
#' Journal of Law and Economics.  
#' 
#' Harris, Curtis C. 1964. “A scientific method of districting”. 
#' Behavioral Science 3(9), 219–225.
#' 
#' Maceachren, A. 1985. Compactness of Geographic Shape: Comparison and 
#' Evaluation of Measures. Geografiska Annaler. Series B, Human Geography, 67(1), 
#' 53-67.
#'
#' Polsby, Daniel D., and Robert D. Popper. 1991. “The Third Criterion: 
#' Compactness as a procedural safeguard against partisan gerrymandering.” 
#' Yale Law & Policy Review 9 (2): 301–353.
#'  
#' Reock, E. 1961. A Note: Measuring Compactness as a Requirement of Legislative 
#' Apportionment. Midwest Journal of Political Science, 5(1), 70-74.
#' 
#' Schwartzberg, Joseph E. 1966. Reapportionment, Gerrymanders, and the Notion 
#' of Compactness. Minnesota Law Review. 1701.
#' 
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom sf st_cast st_bbox st_centroid st_within st_point_on_surface st_coordinates
#' @importFrom sf st_linestring st_intersection st_area st_crs st_is_longlat st_length
#' @importFrom sf st_convex_hull st_crs<- st_geometry st_distance st_union st_touches
#' @importFrom lwgeom st_perimeter st_minimum_bounding_circle
#' @importFrom dplyr select all_of arrange bind_rows rename summarize
#' @importFrom stats dist
#' 
#' @examples
#' \dontrun{
#' data("fl25")
#' data("algdat.p10")
#' 
#' redist.compactness(shp = fl25, district_membership = algdat.p10$cdmat[,1:3], 
#' measure = c('PolsbyPopper', 'EdgesRemoved'))
#'
#' }
#' 
#' @export redist.compactness
redist.compactness <- function(shp = NULL, 
                               district_membership, 
                               measure = c("PolsbyPopper"),
                               population = NULL, adjacency = NULL, nloop = 1,
                               ncores = 1, counties = NULL, pprcpp){
  
  # Check Inputs
  if(is.null(shp)&is.null(adjacency)){
    stop('Please provide a shp or adjacency argument.')
  }  
  
  if(!is.null(shp)){  
    if('SpatialPolygonsDataFrame' %in% class(shp)){
      shp <- shp %>%  st_as_sf()
    } else if(!('sf' %in% class(shp))){
      stop('Please provide "shp" as a SpatialPolygonsDataFrame or sf object.')
    }
  }
  
  if(any(class(district_membership) == 'redist')){
    district_membership <- district_membership$partitions
  }  
  
  if(!any(class(district_membership) %in% c('numeric', 'integer', 'matrix'))){
    stop('Please provide "district_membership" as a numeric vector or matrix.')
  }
  


  
  if("all" %in% measure){
    measure <-  c("PolsbyPopper", "Schwartzberg", "LengthWidth", "ConvexHull", "Reock", "BoyceClark", "FryerHolden", "EdgesRemoved", "logSpanningTree")
  }
  match.arg(arg = measure,several.ok = TRUE, choices = c("PolsbyPopper", "Schwartzberg", "LengthWidth", "ConvexHull", "Reock", "BoyceClark", "FryerHolden", "EdgesRemoved", "logSpanningTree"))
  
  if('FryerHolden' %in% measure & is.null(population)) {
    stop('Please provide a "population" argument when FryerHolden is specified.')
  }  
  
  if('FryerHolden' %in% measure){
    if(!any(class(population) %in% c('numeric', 'integer'))) {
      stop('Please provide "population" as a numeric or integer.')
    }}
  
  if(class(nloop) != 'numeric'){
    stop('Please provide "nloop" as a numeric.')
  }
  
  if(class(ncores) != 'numeric'){
    stop('Please provide "ncores" as a numeric.')
  }
  if(('logSpanningTree' %in% measure) & is.null(counties)){
    stop('Please provide "counties" as an argument.')
  }
  
  
  # Compute compactness scores
  dists <- sort(unique(c(district_membership)))
  nd <-  length(dists)
  
  
  if(!is.matrix(district_membership)){
    district_membership <- as.matrix(district_membership)
  }  
  
  if(missing(pprcpp)){
    if(ncol(district_membership) > 8){
      pprcpp <- TRUE
    } else{
      pprcpp <- FALSE
    }
  }
  
  nmap <-  ncol(district_membership)
  if(nmap!=1){
    nloop = rep(nloop + (1:ncol(district_membership)) - 1, each = nd)
  } else {
    nloop = rep(nloop, nd)
  }
  
  # Initialize object
  comp <- tibble(districts = rep(x = dists, nmap), 
                 PolsbyPopper = rep(NA_real_, nd*nmap), 
                 Schwartzberg = rep(NA_real_, nd*nmap),
                 LengthWidth = rep(NA_real_, nd*nmap),  
                 ConvexHull = rep(NA_real_, nd*nmap),
                 Reock = rep(NA_real_,nd*nmap), 
                 BoyceClark = rep(NA_real_, nd*nmap), 
                 FryerHolden = rep(NA_real_, nd*nmap),
                 EdgesRemoved = rep(NA_real_, nd*nmap), 
                 logSpanningTree = rep(NA_real_, nd*nmap),
                 nloop = nloop) %>% 
    dplyr::select(all_of(c("districts", measure)), all_of(measure), nloop)
  
  if(any(measure %in% c("PolsbyPopper", "Schwartzberg", "LengthWidth", 
                        "ConvexHull", "Reock", "BoyceClark"))){
    nc <- min(ncores, ncol(district_membership))
    if (nc == 1){
      `%oper%` <- `%do%`
    } else {
      `%oper%` <- `%dopar%`
      cl <- makeCluster(nc, setup_strategy = 'sequential')
      registerDoParallel(cl)
      on.exit(stopCluster(cl))
    }
  }
  
  if('PolsbyPopper' %in% measure & pprcpp){
    if(missing(adjacency)){
      warning("For PolsbyPopper, it's recommended that you supply an adjacency list.")
      alist <- lapply(redist.adjacency(shp), function(x){x+1L})
    } else{
      alist <- lapply(adjacency, function(x){x+1L})
    }
    
    invalid <- which(!st_is_valid(shp))
    st_geometry(shp[invalid,]) <- st_geometry(st_buffer(shp[invalid,],0))
    suppressWarnings(perims <- st_perimeter(shp))
    
    perim_adj_df <- lapply(1:length(alist), function(from) {
      suppressWarnings(lines <- st_intersection(shp[from,], shp[alist[[from]],]))
      #suppressWarnings(lines <- st_cast(lines))
      l_lines <- st_length(lines)
      data.frame(origin = from,
                 touching = alist[[from]],
                 edge = as.vector(l_lines))
    }) %>% 
      bind_rows() %>% 
      filter(edge > 0)
    
    adj_boundary_lengths <- perim_adj_df %>% 
      group_by(origin) %>% 
      summarize(perim_adj = sum(edge)) %>%
      mutate(perim_full = as.numeric(perims),
             perim_boundary = perim_full - perim_adj,
             X1 = -1) %>%
      filter(perim_boundary > .001) %>%
      dplyr::select(X1, origin, perim_boundary) %>%
      rename(origin = X1, touching = origin, edge = perim_boundary)
    perim_df <- bind_rows(perim_adj_df, adj_boundary_lengths) %>%
      arrange(origin, touching)
    
    splits <- split(x = district_membership, 
                    rep(1:nc, 
                        each = ceiling(ncol(district_membership)/nc)*nrow(district_membership))[1:(ncol(district_membership)*nrow(district_membership))]) %>% 
      lapply(., FUN = function(x, r = nrow(district_membership)){matrix(data = x, nrow = r)})
    
    result <- foreach(map = 1:nc, .combine = 'cbind', .packages = c('sf', 'lwgeom')) %oper% {
      polsbypopper(from = perim_df$origin, to = perim_df$touching, area = st_area(shp),
                   perimeter = perim_df$edge, dm = splits[[map]], nd = nd)
    }
    
    comp[['PolsbyPopper']] <- c(result)
    
  }
  
  # Compute Specified Scores for provided districts
  if(any(measure %in% c("Schwartzberg", "LengthWidth", 
                        "ConvexHull", "Reock", "BoyceClark") | ('PolsbyPopper' %in% measure & !pprcpp)) ){
    
    nm <- sum(measure %in% c("Schwartzberg", "LengthWidth",
                             "ConvexHull", "Reock", "BoyceClark"))
    
    if('PolsbyPopper' %in% measure & !pprcpp){
      nm <- nm + 1
    }
    
    
    results <- foreach(map = 1:nmap, .combine = 'rbind', .packages = c('sf', 'lwgeom')) %oper% {
      ret <- matrix(nrow = nd, ncol = nm)
      for (i in 1:nd){
        col <- 1
        united <- st_union(shp[district_membership[, map] == dists[i],])
        area <- st_area(united)
        
        if(is.null(st_crs(united$EPSG))||is.na(st_is_longlat(united))){
          perim <- sum(st_length(st_cast(st_cast(united, 'POLYGON'),'LINESTRING')))
        } else if (st_is_longlat(united)){
          perim <- sum(st_length(united))
        } else {
          perim <- sum(st_perimeter(united))
        }
        
        if('PolsbyPopper' %in% measure & !pprcpp){
          ret[i, col] <- 4*pi*(area)/(perim)^2
          col <- col + 1
        }
        
        if('Schwartzberg' %in% measure){
          ret[i, col] <- (perim/(2*pi*sqrt(area/pi)))
          col <- col + 1
        }
        
        if('LengthWidth' %in% measure){
          bbox <- st_bbox(united)
          ratio <- unname((bbox$xmax - bbox$xmin)/(bbox$ymax - bbox$ymin))
          ret[i, col] <- ifelse(ratio>1, 1/ratio, ratio)
          col <- col + 1
        }
        
        if('ConvexHull' %in% measure){
          cvh <- st_area(st_convex_hull(united))
          ret[i, col] <- area/cvh
          col <- col + 1
        }
        
        if('Reock' %in% measure){
          mbc <- st_area(st_minimum_bounding_circle(united))
          ret[i, col] <- area/mbc
          col <- col + 1
        }
        
        if('BoyceClark' %in% measure){
          suppressWarnings(center <- st_centroid(united))
          suppressWarnings(if(!st_within(united,center,sparse=F)[[1]]){
            center <- st_point_on_surface(united)
          })
          center <- st_coordinates(center)
          bbox <- st_bbox(united)
          max_dist <- sqrt((bbox$ymax-bbox$ymin)^2+(bbox$xmax-bbox$xmin)^2)
          st_crs(united) <- NA
          
          x_list <- center[1] + max_dist*cos(seq(0,15)*pi/8)
          y_list <- center[2] + max_dist*sin(seq(0,15)*pi/8)
          radials <- rep(NA_real_, 16)
          for(angle in 1:16){
            line <- data.frame(x = c(x_list[angle],center[1]), y = c(y_list[angle], center[2])) %>%
              st_as_sf(coords = c('x','y'))  %>% st_coordinates() %>% st_linestring()
            radials[angle] <- max(0, dist(st_intersection(line, united)))
          }
          ret[i, col] <- 1 - (sum(abs(radials/sum(radials)*100-6.25))/200)
          col <- col + 1
        }
      }
      return(ret)
    }
    if('PolsbyPopper' %in% measure & pprcpp){
      comp[,3:(nm+2)] <- results
    } else {
      comp[,2:(nm+1)] <- results
    }
  }
  
  if('FryerHolden' %in% measure){
    suppressWarnings(centroids <- st_geometry(st_centroid(shp)))
    dist_sqr <- st_distance(centroids, centroids)^2
    pop <- population*t(matrix(rep(population,nrow(shp)),nrow(shp)))
    fh <- pop*dist_sqr
    comp['FryerHolden'] <- rep(unlist(lapply(1:ncol(district_membership), function(x){
      sum <- sum(unlist(
        lapply(1:nd, function(i){
          ind <- district_membership[,x] == i
          return(sum(fh[ind,ind]))
        })
      ))
      return(sum)
    })), each = nd)
    
  }
  
  
  
  
  if(('EdgesRemoved' %in% measure ||'logSpanningTree' %in% measure) & is.null(adjacency)){
    adjacency <- redist.adjacency(shp)
  }
  
  if('logSpanningTree' %in% measure){
    comp[['logSpanningTree']] <- rep(log_st_map(g = adjacency, 
                                                districts = district_membership, 
                                                counties = counties,
                                                n_distr = nd), each = nd)
  }
  
  if('EdgesRemoved' %in% measure){
    comp[['EdgesRemoved']] <- rep(n_removed(g = adjacency, 
                                            districts = district_membership, 
                                            n_distr = nd), each = nd)
  }
  
  # Return results
  return(comp)
}

utils::globalVariables(c("i", "j", "edge", "origin", "perim_full", "perim_adj",
                         "perim_boundary", "X1", ".", "touching"))
