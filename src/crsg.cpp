#include <RcppArmadillo.h>
#include "redist_types.h"
#include "distance_helpers.h"
#include "make_swaps_helper.h"
#include "constraint_calc_helper.h"
#include "shatter_search.h"
using namespace Rcpp;

// [[Rcpp::export]]
List crsg(List adj_list,
          NumericVector population,
          NumericVector area,
          NumericVector x_center,
          NumericVector y_center,
          int Ndistrict,
          double target_pop,
          double thresh,
          int maxiter){
  
  // Initializations
  int n_prec = adj_list.size();
  int i_dist, j_dist, temp_dist, i = 0, j = 0, prec_swap;
  NumericVector dist;
  Rcpp::List adj_list_up = clone(adj_list), adj_list_temp = clone(adj_list);
  NumericVector x = clone(x_center), y = clone(y_center), area_up = clone(area), district_pop(Ndistrict);
  List member_plist(n_prec);
  IntegerVector prec(n_prec), mp_i, mp_j, adj_i, adj_j;
  IntegerVector empty_vec = IntegerVector::create(NA_INTEGER);
  bool pop_swap, no_ij, to_next, cont_swap;
  int size = 0;
  int HPP, LPP;
  NumericVector new_dist_assignment(n_prec);
  List alConnected;
  
  
  // Step 1
  int n_current = n_prec;
  IntegerVector dist_assignment = seq(0, n_prec - 1);
  for(i = 0; i < n_prec; i++){
    member_plist[i] = dist_assignment[i];
  }
  IntegerVector dist_sample = dist_assignment;
  
  
  // Step 2
  while(n_current > Ndistrict){
    
    // Step 2A:
    i_dist = as<int>(sample(dist_sample, 1));
    
    // Step 2B:
    adj_i = adj_list_up[i_dist];
    j_dist = closest_adj(adj_i, i_dist, x, y);
    
    // Step 2C:
    // ensure i < j for consistency
    if(i_dist > j_dist){
      temp_dist = i_dist;
      i_dist = j_dist;
      j_dist = temp_dist;
    }
    
    // update i dist information
    // 1 - merge the member lists of i and j
    mp_i = member_plist[i_dist];
    mp_j = member_plist[j_dist];
    member_plist[i_dist] = union_(mp_i, mp_j);
    // 2 - merge the adjacencies of i and j
    adj_i = adj_list_up[i_dist];
    adj_j = adj_list_up[j_dist];
    adj_i = union_(adj_i, adj_j);
    adj_list_up[i_dist] = unique(adj_i);
    // 3.1 - Calculate new center of gravity by weighting
    x[i_dist] = (x[i_dist]*area_up[i_dist] + x[j_dist]*area_up[j_dist])/(area_up[i_dist]+area_up[j_dist]);
    y[i_dist] = (y[i_dist]*area_up[i_dist] + y[j_dist]*area_up[j_dist])/(area_up[i_dist]+area_up[j_dist]);
    // 3.2 - update the new area weights
    area_up[i_dist] = area_up[i_dist] + area_up[j_dist];
    
    // remove j dist information
    // 1  - remove j plist
    member_plist.erase(j_dist);
    // 2 - erase j adjacency
    adj_list_up.erase(j_dist);
    // 3 - erase the j area
    area_up.erase(j_dist);
    // 4 - erase the centers 
    x.erase(j_dist);
    y.erase(j_dist);
    
    // Sink Adjacency
    for(i = 0; i < adj_list_up.size(); i++){
      adj_i = adj_list_up[i];
      for(j = 0; j < adj_i.size(); j++){
        if(adj_i[j] == j_dist){
          adj_i[j] = i_dist;
        }else if(adj_i[j] > j_dist){
          adj_i[j] = adj_i[j] - 1;
        }
      }
      adj_list_up[i] = unique(adj_i);
    }
    // Remove self adjacency
    adj_i = adj_list_up[i_dist];
    j = adj_i.size(); i = 0;
    while(j == adj_i.size() && i < adj_i.size()){
      if(adj_i[i] == i_dist){
        j = i;
      }
      i++;
    }
    adj_i.erase(j);
    adj_list_up[i_dist] = adj_i;
    
    // Update district assignment
    for(i = 0; i < n_prec; i++){
      if(dist_assignment[i] == j_dist){
        dist_assignment[i] = i_dist;
      }
      if(dist_assignment[i] > j_dist){
        dist_assignment[i] = dist_assignment[i] - 1;
      }
    }
    
    // update loop params
    n_current = n_current - 1;
    dist_sample.erase(n_current);
  }
  
  // Create parameters for step 3
  // create population list
  for(i = 0; i < n_prec; i++){
    district_pop[dist_assignment[i]] += population[i];
  }
  // create population bounds
  double min_pop = target_pop*(1-thresh);
  double max_pop = target_pop*(1+thresh);
  // create exit condition
  int iter = 0;

  // Step 3
  while((max(district_pop) > max_pop) || (min(district_pop) < min_pop) ){
    
    // First build the right size vector
    size = 0;
    for(i = 0; i < Ndistrict; i++){
      adj_i = adj_list_up[i];
      size += adj_i.size();
    }
    
    // Create a population by district difference set of vectors
    IntegerVector pop_dist(size), i_vec(size), j_vec(size);
    size = 0;
    for(i = 0; i < Ndistrict; i++){
      adj_i = adj_list_up[i];
      for(j = 0; j < adj_i.size(); j++){
        pop_dist[size] = district_pop[i] - district_pop[ adj_i[j] ];
        i_vec[size] = i;
        j_vec[size] = adj_i[j];
        size++;
      }
    }
    
    // Rep code only allows swap if it doesnt make pop parity worse
    LPP = min(district_pop);
    HPP = max(district_pop);
    
    // Rep code iterates over remaining ij districts if no swap available in optimal district
    no_ij = TRUE;
    while(no_ij){
      
      // Rare condition where no swaps are possible without making things worse
      // implies threshold too strong for this random seed
      if(pop_dist.size() == 0){
        List result;
        result["plan"] = NumericVector::get_na();
        result["district_list"] = NumericVector::get_na();
        result["district_pop"] = NumericVector::get_na();
        return(result);
        }
      
      // Pick the vector with the worse population difference
      i_dist = i_vec[which_max(pop_dist)];
      j_dist = j_vec[which_max(pop_dist)];
      

      // Step 3B
      mp_i = member_plist[i_dist];
      IntegerVector swappable(mp_i.size());
      for(i = 0; i < mp_i.size(); i++){
        adj_i = adj_list[mp_i[i]];
        if(can_swap(adj_list, mp_i[i], j_dist, dist_assignment)){
          swappable[i] = 1;
        }
      }
      
      IntegerVector swap_candidates(sum(swappable));
      j = 0;
      for(i = 0; i < mp_i.size(); i++){
        if(swappable[i] == 1){
          swap_candidates[j] = mp_i[i];
          j++;
        }
      }

      // Step 3C
      // Rep code calculates all of these to function as a ranking in 3D
      NumericVector dist(swap_candidates.size());
      size = swap_candidates.size();
      for(i = 0; i < swap_candidates.size(); i++){
        dist[i] = dist_dist_diff(swap_candidates[i], i_dist, j_dist,
                                 x_center, y_center, x, y); 
      }

      // Step 3D
      // if there is at least one possible entry, we continue
      if(dist.size() > 0){
        prec_swap = swap_candidates[which_max(dist)];
        
        // Step 3D.5 comes from rep code too
        // we loop over the other precincts in the district if the first option makes
        // population parity worse; note that should then avoid infinite loops
        pop_swap = ((district_pop[j_dist] + population[prec_swap]) < HPP) &&
          ((district_pop[i_dist] - population[prec_swap]) > LPP );
        
        
        // Check contiguity here:
        new_dist_assignment = as<NumericVector>(clone(dist_assignment));
        new_dist_assignment(prec_swap) = (double)j_dist;
        alConnected = genAlConn(adj_list, new_dist_assignment);
        cont_swap = countpartitions(alConnected) == Ndistrict;
        
        to_next = TRUE;
        while(!(pop_swap && cont_swap) && to_next){
          // If there are at least two choices proceed as normal
          if(dist.size() > 1){
            swap_candidates.erase(which_max(dist));
            dist.erase(which_max(dist));
            prec_swap = swap_candidates[which_max(dist)];
            pop_swap = ((district_pop[j_dist] + population[prec_swap]) < HPP) &&
              ((district_pop[i_dist] - population[prec_swap]) > LPP );
            
            if(pop_swap){
              new_dist_assignment = as<NumericVector>(clone(dist_assignment));
              new_dist_assignment(prec_swap) = (double)j_dist;
              alConnected = genAlConn(adj_list, new_dist_assignment);
              cont_swap = countpartitions(alConnected) == Ndistrict;
            }
            
            
            if(pop_swap && cont_swap){//exp
              no_ij = FALSE;
              to_next = FALSE;
            }
          } else{ // then only one prec remains for this ij
            // Try it
            prec_swap = swap_candidates[0];
            pop_swap = ((district_pop[j_dist] + population[prec_swap]) < HPP) &&
              ((district_pop[i_dist] - population[prec_swap]) > LPP );
            to_next = FALSE;
            // If it won't work move to next pair ij
            if(pop_swap){
              new_dist_assignment = as<NumericVector>(clone(dist_assignment));
              new_dist_assignment(prec_swap) = (double)j_dist;
              alConnected = genAlConn(adj_list, new_dist_assignment);
              cont_swap = countpartitions(alConnected) == Ndistrict;
            }
            
            to_next = FALSE;
            // If it won't work move to next pair ij
            if(!( pop_swap && cont_swap  )){
              no_ij = TRUE;
              i_vec.erase(which_max(pop_dist));
              j_vec.erase(which_max(pop_dist));
              pop_dist.erase(which_max(pop_dist));
            } else{ // if it works just update the loop logic
              no_ij = FALSE;
            }
          }
        }
        
        // if we have found something and made it here, we want to update the loop
        if(pop_swap){
          no_ij = FALSE;
        }
        
      } else {// end if dist.size() > 0
        // if dist size is zero, we need next ij
        no_ij = TRUE;
        i_vec.erase(which_max(pop_dist));
        j_vec.erase(which_max(pop_dist));
        pop_dist.erase(which_max(pop_dist));
      }
    } // end while(no ij) 

    // perform the swap
    // add to member list
    mp_j = member_plist[j_dist];
    mp_j.push_back(prec_swap);
    member_plist[j_dist] = mp_j;
    // remove from old list
    mp_i = member_plist[i_dist];
    for(i = 0; i < mp_i.size(); i++){
      if(mp_i[i] == prec_swap){
        mp_i.erase(i);
      }
    }
    member_plist[i_dist] = mp_i;
    
    // update in dist assignment
    dist_assignment[prec_swap] = j_dist;
    
    // Update corresponding variables
    // Adjust centers of gravity:
    x[i_dist] = (x[i_dist]*area_up[i_dist] - x_center[prec_swap]*area[prec_swap])/(area_up[i_dist] - area[prec_swap]);
    y[i_dist] = (y[i_dist]*area_up[i_dist] - y_center[prec_swap]*area[prec_swap])/(area_up[i_dist] - area[prec_swap]);
    x[j_dist] = (x[j_dist]*area_up[j_dist] + x_center[prec_swap]*area[prec_swap])/(area_up[j_dist] + area[prec_swap]);
    y[j_dist] = (y[j_dist]*area_up[j_dist] + y_center[prec_swap]*area[prec_swap])/(area_up[j_dist] + area[prec_swap]);
    // Adjust Areas:
    area_up[i_dist] -= area[prec_swap];
    area_up[j_dist] += area[prec_swap];
    // Adjust populations
    district_pop[i_dist] -= population[prec_swap];
    district_pop[j_dist] += population[prec_swap];
    // Recalculate district adjacencies:
    adj_list_temp = clone(adj_list);
    for(i = 0; i < adj_list_temp.size(); i++){
      adj_i = adj_list_temp[i];
      for(j = 0; j < adj_i.size(); j++){
        adj_i[j] = dist_assignment[adj_i[j]];
      }
      adj_list_temp[i] = unique(adj_i);
    }
    for(i = 0; i < Ndistrict; i++){
      adj_list_up[i] = empty_vec;
    }
    for(i = 0; i < n_prec; i++){
      i_dist = dist_assignment[i];
      adj_i = adj_list_up[i_dist];
      adj_j = adj_list_temp[i];
      adj_list_up[i_dist] = union_(adj_i, adj_j);
    }
    for(i = 0; i < Ndistrict; i++){
      adj_i = adj_list_up[i];
      adj_i = na_omit(adj_i);
      for(j = 0; j < adj_i.size(); j++){
        if(adj_i[j] == i){
          adj_i.erase(j);
        }
      }
      adj_list_up[i] = adj_i;
      
    }
    
    // End if program has gone too long to avoid endless loops
    iter++;
    if(iter > maxiter){
      List result;
      result["plan"] = NumericVector::get_na();
      result["district_list"] = NumericVector::get_na();
      result["district_pop"] = NumericVector::get_na();
      return(result);
    }
  }
  
  
  // Output is:
  // member_dvec - Vector of length #precincts with CD assigments
  // member_plist - List of length #districts with precinct ID's
  // district_pop - vector of length #districts with district populations
  List result;
  result["plan"] = dist_assignment;
  result["district_list"] = member_plist;
  result["district_pop"] = district_pop;
  return(result);
  // Output structure directly mimics RSG output from earlier releases
  
}
