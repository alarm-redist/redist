///////////////////////////////////////////////
// Author: Ben Fifield
// Institution: Princeton University
// Date Created: 2014/12/26
// Date Last Modified: 2015/02/26
// Purpose: Contains functions to run calculate beta constraints
/////////////////////////////////////////////// 

// Header files
#include <RcppArmadillo.h>

using namespace Rcpp;

/* Function to identify which precincts lie on the boundary of a congressional
   district */
// [[Rcpp::export]]
NumericVector findBoundary(List fullList,
			   List conList)
{

  /* Inputs to function:
     fullList: Full adjacency list of geographic units

     conList: Adjacency list of geographic units within cong district
  */

  // Initialize container vector of 0's (not boundary) and 1's (boundary)
  NumericVector isBoundary(fullList.size());

  // Initialize inside loop
  NumericVector full; NumericVector conn; int i;

  // Loop through aList
  for(i = 0; i < fullList.size(); i++){

    // Get vectors of full and cd-connected components for precinct i
    full = fullList(i);
    conn = conList(i);

    // Compare lengths - if conn < full, then boundary unit
    if(full.size() > conn.size()){
      isBoundary(i) = 1;
    }
    
  }

  return isBoundary;

}

// Fryer-Holden measure
List fh_compact(arma::uvec new_cds,
		arma::uvec current_cds,
		NumericVector pops,
		NumericMatrix ssdmat,
		double denominator = 1.0){

  // Initialize objects
  double ssd_new = 0.0;
  double ssd_old = 0.0;
  int j; int k;

  // SSD for new partition
  for(j = 0; j < new_cds.size(); j++){
    for(k = j + 1; k < new_cds.size(); k++){
      ssd_new += (double)ssdmat(new_cds(j),new_cds(k)) *
	pops(new_cds(j)) * pops(new_cds(k));
    }
  }

  // SSD for old partition
  for(j = 0; j < current_cds.size(); j++){
    for(k = j + 1; k < current_cds.size(); k++){
      ssd_old += (double)ssdmat(current_cds(j),current_cds(k)) *
	pops(current_cds(j)) * pops(current_cds(k));
    }
  }

  List out;
  out["ssd_new"] = ssd_new / denominator;
  out["ssd_old"] = ssd_old / denominator;
  
}

// Polsby-popper measure
List pp_compact(arma::uvec new_cds,
		arma::uvec current_cds,
		NumericVector areas_vec,
		NumericVector boundarylist_new,
		NumericVector boundarylist_current,
		List aList,
		List boundarylength_list){

  // Declare objects
  arma::uvec new_boundaryprecs = find(boundarylist_new == 1);
  arma::uvec current_boundaryprecs = find(boundarylist_current == 1);
  arma::ivec new_boundaryprecs_indist = intersect(new_cds, new_boundaryprecs);
  arma::ivec current_boundaryprecs_indist = intersect(current_cds, current_boundaryprecs);
  
  double area_new = 0.0;
  double area_old = 0.0;
  double perimeter_new = 0.0;
  double perimeter_old = 0.0;
  int j; int k;

  arma::vec perimeter_vec;
  arma::vec adj_precs;
  arma::vec adj_boundary;
  arma::uvec adj_precs_gt;
  arma::uvec adj_precs_inds;
  arma::uvec new_boundaryprecs_indist_inds;
  arma::ivec indices_boundary_indist;

  double pi = 3.141592653589793238463;

  // Areas for current and new partitions
  for(j = 0; j < new_cds.size(); j++){
    area_new += areas_vec(new_cds(j));
  }
  for(j = 0; j < current_cds.size(); j++){
    area_old += areas_vec(current_cds(j));
  }

  // Perimeters for current and new partitions
  for(j = 0; j < new_boundaryprecs_indist.n_elem; j++){
    adj_precs = aList(new_boundaryprecs_indist[j]);
    perimeter_vec = boundarylength_list(newboundaryprecs_indist[j]);
    // Get indices of adj_precs that are greater than j
    // adj_precs_gt
    adj_precs_gt = find(adj_precs > j);
    // Get indices of adj_precs that are boundaries
    // adj_precs_inds
    intersect(adj_boundary, adj_precs_inds, new_boundaryprecs_indist_inds, adj_precs, new_boundaryprecs_indist);
    // Get intersection
    indices_boundary_indist = intersect(adj_precs_gt, adj_precs_inds);
    for(k = 0; k < indices_boundary_indist.n_elem; k++){
      perimeter_new += perimeter_vec(indices_boundary_indist[k]);
    }
  }
  for(j = 0; j < current_boundaryprecs_indist.n_elem; j++){
    adj_precs = aList(current_boundaryprecs_indist[j]);
    perimeter_vec = boundarylength_list(currentboundaryprecs_indist[j]);
    // Get indices of adj_precs that are greater than j
    // adj_precs_gt
    adj_precs_gt = find(adj_precs > j);
    // Get indices of adj_precs that are boundaries
    // adj_precs_inds
    intersect(adj_boundary, adj_precs_inds, current_boundaryprecs_indist_inds, adj_precs, current_boundaryprecs_indist);
    // Get intersection
    indices_boundary_indist = intersect(adj_precs_gt, adj_precs_inds);
    for(k = 0; k < indices_boundary_indist.n_elem; k++){
      perimeter_old += perimeter_vec(indices_boundary_indist[k]);
    }
  }

  List out;
  out["pp_new"] = (double)4.0 * pi * area_new / pow(perimeter_new, 2.0);
  out["pp_old"] = (double)4.0 * pi * area_old / pow(perimeter_old, 2.0);
  
}

// Function to calculate the strength of the beta constraint for population
List calc_psipop(arma::vec current_dists,
		 arma::vec new_dists,
		 NumericVector pops,
		 NumericVector distswitch)
{

  /* Inputs to function 
     current_dists: vector of the current cong district assignments
     new_dists: vector of the new cong district assignments
     pops: vector of district populations
     weight_population: strength of the beta constraint
     distswitch: vector containing the old district, and the proposed new district
  */

  // Calculate parity
  double parity = (double)sum(pops) / (max(current_dists) + 1);

  // Initialize psi values
  double psi_new = 0.0;
  double psi_old = 0.0;

  // Loop over congressional districts
  for(int i = 0; i < distswitch.size(); i++){

    // Population objects
    int pop_new = 0;
    int pop_old = 0;
    arma::uvec new_cds = find(new_dists == distswitch(i));
    arma::uvec current_cds = find(current_dists == distswitch(i));

    // Get population of the old districts
    for(int j = 0; j < new_cds.size(); j++){
      pop_new += pops(new_cds(j));
    }
    for(int j = 0; j < current_cds.size(); j++){
      pop_old += pops(current_cds(j));
    }

    // Calculate the penalty
    psi_new += std::pow(pop_new / parity - 1.0, 2.0);
    psi_old += std::pow(pop_old / parity - 1.0, 2.0);

  }

  // Create return object
  List out;
  out["pop_new_psi"] = std::sqrt(psi_new);
  out["pop_old_psi"] = std::sqrt(psi_old);

  return out;

}

// Function to calculate the strength of the beta constraint for compactness
// Currently implemented: Fryer and Holden 2011 RPI index, Polsby-Popper
List calc_psicompact(arma::vec current_dists,
		     arma::vec new_dists,
		     NumericVector distswitch,
		     std::string measure,
		     // For Polsby-Popper
		     List aList,
		     NumericVector areas_vec,
		     List boundarylength_list,
		     // For Fryer Holden
		     NumericVector pops,
		     NumericMatrix ssdmat,
		     double denominator = 1.0){

  /* Inputs to function:
     current_dists: vector of the current cong district assignments
     new_dists: vector of the new cong district assignments
     pops: vector of district populations
     beta_compact: strength of the beta constraint
     distswitch: vector containing the old district, and the proposed new district
     ssdmat: squared distance matrix
     denominator: normalizing constant for rpi
  */
  
  // Initialize psi values
  double psi_new = 0.0;
  double psi_old = 0.0;

  // Initialize lists and boundary vectors
  List aList_new;
  List aList_current;
  NumericVector boundarylist_new(new_dists.size());
  NumericVector boundarylist_current(current_dists.size());
  if(measure == "polsby-popper"){
    aList_new = genAlConn(aList, new_dists);
    aList_current = genAlConn(aList, current_dists);
    boundarylist_new = findBoundary(aList, aList_new);
    boundarylist_current = findBoundary(aList, aList_current);
  }

  // Loop over the congressional districts
  for(int i = 0; i < distswitch.size(); i++){

    // Initialize objects
    arma::uvec new_cds = find(new_dists == distswitch(i));
    arma::uvec current_cds = find(current_dists == distswitch(i));

    if(measure == "fryer-holden"){

      List fh_out = fh_compact(new_cds, current_cds, pops, ssdmat, denominator);
      
      // Add to psi
      psi_new += fh_out["ssd_new"];
      psi_old += fh_out["ssd_old"];
      
    }else if(measure == "polsby-popper"){

      List pp_out = pp_compact(new_cds, current_cds, areas_vec, boundarylist_new,
			       boundarylist_current, aList, boundarylength_list);

	// Add to psi
      psi_new += pp_out["pp_new"];
      psi_old += pp_out["pp_old"];
      
    }

  }

  // Create return object
  List out;
  out["compact_new_psi"] = psi_new;
  out["compact_old_psi"] = psi_old;

  return out;

}

// Function to constrain by segregating a group
List calc_psisegregation(arma::vec current_dists,
			 arma::vec new_dists,
			 NumericVector pops,
			 NumericVector distswitch,
			 NumericVector grouppop)
{

  /* Inputs to function:
     current_dists: vector of the current cong district assignments
     new_dists: vector of the new cong district assignments
     pops: vector of district populations
     beta_segregation: strength of the beta constraint
     distswitch: vector containing the old district, and the proposed new district
     grouppop: vector of subgroup district populations
     
  */

  // Initialize psi values
  double psi_new = 0.0;
  double psi_old = 0.0;

  // Initialize denominator
  int T = sum(pops);
  double pAll = (double)sum(grouppop) / T;
  double denom = (double)2 * T * pAll * (1 - pAll);
  
  // Loop over congressional districts
  for(int i = 0; i < distswitch.size(); i++){

    // Initialize objects
    int oldpopall = 0;
    int newpopall = 0;
    int oldpopgroup = 0;
    int newpopgroup = 0;
    arma::uvec new_cds = find(new_dists == distswitch(i));
    arma::uvec current_cds = find(current_dists == distswitch(i));
  
    // Segregation for proposed assignments
    for(int j = 0; j < new_cds.size(); j++){
      newpopall += pops(new_cds(j));
      newpopgroup += grouppop(new_cds(j));
    }
  
    // Segregation for current assignments
    for(int j = 0; j < current_cds.size(); j++){
      oldpopall += pops(current_cds(j));
      oldpopgroup += grouppop(current_cds(j));
    }
  
    // Calculate proportions
    // Rcout << "old population group " << oldpopgroup << std::endl;
    // Rcout << "old population all " << oldpopall << std::endl;
    double oldgroupprop = (double)oldpopgroup / oldpopall;
    // Rcout << "old proportion group " << oldgroupprop << std::endl;
    double newgroupprop = (double)newpopgroup / newpopall;

    // Get dissimilarity index
    psi_new += (double)(newpopall * std::abs(newgroupprop - pAll));
    psi_old += (double)(oldpopall * std::abs(oldgroupprop - pAll));

  }
  
  // Standardize psi
  psi_new = (double)psi_new / denom;
  psi_old = (double)psi_old / denom;

  // Create return object
  List out;
  out["segregation_new_psi"] = psi_new;
  out["segregation_old_psi"] = psi_old;

  return out;

}

// Function to constrain on plan similarity to original plan
List calc_psisimilar(arma::vec current_dists,
		     arma::vec new_dists,
		     arma::vec orig_dists,
		     NumericVector distswitch)
{

  /* Inputs to function:
     current_dists: vector of the current cong district assignments
     new_dists: vector of the new cong district assignments
     orig_dists: vector of the true congressional district assignments
     beta_similar: strength of the beta constraint
     distswitch: vector containing the old district, and the proposed new district
  */

  // Initialize psi values
  double psi_new = 0.0;
  double psi_old = 0.0;

  // Loop over congressional districts
  for(int i = 0; i < distswitch.size(); i++){

    // Initialize objects
    int new_count = 0;
    int old_count = 0;
    NumericVector orig_cds = wrap(find(orig_dists == distswitch(i)));
    arma::uvec new_cds = find(new_dists == distswitch(i));
    arma::uvec current_cds = find(current_dists == distswitch(i));

    // Similarity measure for proposed assignments
    for(int j = 0; j < new_cds.size(); j++){
      if(any(orig_cds == new_cds(j)).is_true()){
    	new_count++;
      }
    }

    // Similarity measure for current assignments
    for(int j = 0; j < current_cds.size(); j++){
      if(any(orig_cds == current_cds(j)).is_true()){
    	old_count++;
      }
    }

    // Calculate proportions
    double old_count_prop = (double)old_count / orig_cds.size();
    double new_count_prop = (double)new_count / orig_cds.size();
    
    // Add to psi
    psi_new += (double)std::abs(new_count_prop - 1);
    psi_old += (double)std::abs(old_count_prop - 1);

  }

  // Normalize by dividing by number of congressional districts
  psi_new = psi_new / distswitch.size();
  psi_old = psi_old / distswitch.size();

  // Create return object
  List out;
  out["similar_new_psi"] = psi_new;
  out["similar_old_psi"] = psi_old;

  return out;
}

// Function to calculate the strength of the county split penalty
List calc_psicounty(arma::vec current_dists,
		    arma::vec new_dists,
		    arma::vec county_assignments)
{

  // Initialize psi values
  double psi_new = 0.0;
  double psi_old = 0.0;

  /* We want to:
     1) Get the unique county labels
     2) Loop through the unique county labels
     3) Increment psi if there are more than 2 unique district assignments for 
     that county
  */

  // Get unique psi labels
  arma::vec unique_county = unique(county_assignments);

  // Loop over county labels
  int i;
  arma::uvec in_county;
  arma::vec current_dists_incounty;
  arma::vec new_dists_incounty;
  arma::vec unique_current;
  arma::vec unique_new;
  for(i = 0; i < unique_county.n_elem; i++){

    // Get indices that match i
    in_county = find(county_assignments == i);
    current_dists_incounty = current_dists.elem(in_county);
    new_dists_incounty = new_dists.elem(in_county);
    
    // Count unique elements in each
    unique_current = unique(current_dists_incounty);
    unique_new = unique(new_dists_incounty);

    // Increment psi
    if(unique_new.n_elem > 1){
      psi_new += 1.0;
    }
    if(unique_current.n_elem > 1){
      psi_old += 1.0;
    }
    
  }

  psi_new = (double)psi_new / unique_county.n_elem;
  psi_old = (double)psi_old / unique_county.n_elem;

  // Create return object
  List out;
  out["countysplit_new_psi"] = psi_new;
  out["countysplit_old_psi"] = psi_old;

  return out;
  
}

