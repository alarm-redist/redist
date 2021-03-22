flip_constraints_helper <- function(map, 
                                   constraint = 'compact', 
                                   constraintweight = 0.6, 
                                   compactness_metric = 'edges-removed', 
                                   counties,
                                   partisan_metric = 'efficiency-gap',
                                   rvote, 
                                   dvote, 
                                   areas, 
                                   borderlength_mat, 
                                   ssdmat, 
                                   ssd_denom = 1,
                                   group_pop, 
                                   tgt_min = 0.55,
                                   tgt_other = 0.25,
                                   minorityprop
                                   ){
  if(missing(map)){
    stop('Please provide a redist_map object to map.')
  }
  
  
  
  constraints_list <- process_flip_constr(constraints = list(), 
                                         group_pop = rep(0, nrow(map)),
                                         counties = rep(1, nrow(map)))
  

  
}