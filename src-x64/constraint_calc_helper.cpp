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

/* Function to modify adjacency list to reflect adjacency only within
   a particular congressional district */
// [[Rcpp::export]]
List genAlConn(List aList,
	       NumericVector cds)
{

  /* Inputs to function:
     aList: adjacency list of geographic units

     cds: vector of congressional district assignments
  */

  // Initialize container list
  List alConnected(cds.size());

  // Initialize
  int i; NumericVector avec; int cd_i; int j;

  // Loop through precincts
  for(i = 0; i < cds.size(); i++){

    // For precinct i, get adjacent precincts
    avec = aList(i);

    // Get precinct i's congressional district
    cd_i = cds(i);

    // Initialize empty vector
    NumericVector avec_cd;

    // Loop through avec to identify which are in same cd
    for(j = 0; j < avec.size(); j++){

      // Check if j'th entry in avec is same cd, add to avec_cd if so
      if(cds(avec(j)) == cd_i){
	avec_cd.push_back(avec(j));
      }

    }

    // Add to alConnected list
    alConnected(i) = avec_cd;

  }

  return alConnected;

}

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

arma::uvec getIn(arma::ivec vec1, arma::ivec vec2){

  int i; int j; bool match; arma::uvec store_in(vec1.n_elem);
  for(i = 0; i < vec1.n_elem; i++){
    match = false;
    for(j = 0; j < vec2.n_elem; j++){
      if(vec1(i) == vec2(j)){
	match = true;
	break;
      }
    }
    store_in(i) = match;
  }

  return store_in;

}

arma::uvec get_in_index(arma::vec vec1, arma::vec vec2){

  int i; int j; bool match; arma::uvec store_in(vec1.n_elem);
  for(i = 0; i < vec1.n_elem; i++){
    match = false;
    for(j = 0; j < vec2.n_elem; j++){
      if(vec1(i) == vec2(j)){
	match = true;
	break;
      }
    }
    store_in(i) = j;
  }

  return store_in;

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

  return out;

}

// Polsby-popper measure
List pp_compact(arma::uvec new_cds,
		arma::uvec current_cds,
		arma::vec areas_vec,
		arma::vec boundarylist_new,
		arma::vec boundarylist_current,
		arma::mat borderlength_mat,
		arma::vec pop_vec,
		List aList,
		bool discrete = false){

  double pi = 3.141592653589793238463;

  /*
     Calculate the area for the district
     by summing over the areas of the cds
  */
  double area_new = 0.0;
  double area_old = 0.0;
  int j; int k; int l;
  if(discrete == true){
    for(j = 0; j < new_cds.n_elem; j++){
      area_new += pop_vec(new_cds(j));
    }
    for(j = 0; j < current_cds.n_elem; j++){
      area_old += pop_vec(current_cds(j));
    }
  }else{
    for(j = 0; j < new_cds.n_elem; j++){
      area_new += areas_vec(new_cds(j));
    }
    for(j = 0; j < current_cds.n_elem; j++){
      area_old += areas_vec(current_cds(j));
    }
  }

  /*
     Calculate the perimeter for the district
     by finding boundaries on the CD
  */

  // Unpack the borderlength matrix
  arma::vec borderlength_col1 = borderlength_mat.col(0);
  arma::vec borderlength_col2 = borderlength_mat.col(1);
  arma::vec borderlength_col3 = borderlength_mat.col(2);

  // Modify the boundary indicator to include -1s
  arma::uvec border_inds = find(borderlength_col1 == -1.0);
  arma::vec get_unique_borderinds = unique(borderlength_col2.elem(border_inds));
  boundarylist_new.elem(arma::conv_to<arma::uvec>::from(get_unique_borderinds)).replace(0, 1);
  boundarylist_current.elem(arma::conv_to<arma::uvec>::from(get_unique_borderinds)).replace(0, 1);

  // Get indices of the obsevations that are both in the district and
  // on the boundary
  arma::uvec boundary_inds_new = find(boundarylist_new == 1);
  arma::ivec boundary_precs_new = arma::intersect(arma::conv_to<arma::ivec>::from(new_cds), arma::conv_to<arma::ivec>::from(boundary_inds_new));
  arma::uvec boundary_inds_current = find(boundarylist_current == 1);
  arma::ivec boundary_precs_current = arma::intersect(arma::conv_to<arma::ivec>::from(current_cds), arma::conv_to<arma::ivec>::from(boundary_inds_current));

  // Loop over the indices in boundary_precs_new and check the adjacent units for
  // whether they too lie on a boundary
  arma::vec adj_precs;
  arma::uvec in_cd;
  arma::vec adj_precs_sub;

  arma::uvec adj_precs_sub_indices;
  arma::uvec find_boundarylength_in_boundarymat;
  arma::ivec loop_over_inds;
  double perimeter_new = 0.0; double perimeter_old = 0.0;
  for(j = 0; j < boundary_precs_new.n_elem; j++){

    if(discrete == false){

      // Get the adjacent indices - on the boundary and
      // greater than boundary_precs[j] to avoid double counting
      adj_precs = as<arma::vec>(aList(boundary_precs_new(j)));

      // Of the adjacent precincts, which are not in the same cong district?
      in_cd = getIn(arma::conv_to<arma::ivec>::from(adj_precs), arma::conv_to<arma::ivec>::from(new_cds));
      adj_precs_sub = adj_precs.elem( find(in_cd == false) );
      adj_precs_sub.resize(adj_precs_sub.size() + 1);
      adj_precs_sub(adj_precs_sub.size() - 1) = -1.0;

      // Which elements of boundary_mat's first column are in adj_precs_sub,
      // and what of second column equal the actual precinct
      find_boundarylength_in_boundarymat = find(borderlength_col2 == boundary_precs_new(j));
      for(k = 0; k < adj_precs_sub.n_elem; k++){

	adj_precs_sub_indices = find(borderlength_col1 == adj_precs_sub(k));
	loop_over_inds = arma::intersect(arma::conv_to<arma::ivec>::from(adj_precs_sub_indices), arma::conv_to<arma::ivec>::from(find_boundarylength_in_boundarymat));
	if(loop_over_inds.n_elem > 0){
	  for(l = 0; l < loop_over_inds.n_elem; l++){
	    perimeter_new += (double)borderlength_col3(loop_over_inds(l));
	  }
	}

      }

    }else{

      // Add the population around the perimeter
      perimeter_new += (double)pop_vec(boundary_precs_new(j));

    }

  }
  for(j = 0; j < boundary_precs_current.n_elem; j++){

    if(discrete == false){

      // Get the adjacent indices - on the boundary and
      // greater than boundary_precs[j] to avoid double counting
      adj_precs = as<arma::vec>(aList(boundary_precs_current(j)));

      // Of the adjacent precincts, which are not in the same congressional district?
      in_cd = getIn(arma::conv_to<arma::ivec>::from(adj_precs), arma::conv_to<arma::ivec>::from(current_cds));
      adj_precs_sub = adj_precs.elem( find(in_cd == false) );
      adj_precs_sub.resize(adj_precs_sub.size() + 1);
      adj_precs_sub(adj_precs_sub.size() - 1) = -1.0;

      // Which elements of boundary_mat's first column are in adj_precs_sub, and what of second column equal the actual precinct
      find_boundarylength_in_boundarymat = find(borderlength_col2 == boundary_precs_current(j));
      for(k = 0; k < adj_precs_sub.n_elem; k++){

	adj_precs_sub_indices = find(borderlength_col1 == adj_precs_sub(k));
	loop_over_inds = arma::intersect(arma::conv_to<arma::ivec>::from(adj_precs_sub_indices), arma::conv_to<arma::ivec>::from(find_boundarylength_in_boundarymat));
	if(loop_over_inds.n_elem > 0){
	  for(l = 0; l < loop_over_inds.n_elem; l++){
	    perimeter_old += (double)borderlength_col3(loop_over_inds(l));
	  }
	}

      }

    }else{

      perimeter_old += (double)pop_vec(boundary_precs_current(j));

    }

  }

  // Note - multiplying by -1 since we want to maximize Polsby-Popper,
  // but in general minimize metrics
  List out;
  if(discrete == false){
    out["pp_new"] = (double)-1.0 * 4.0 * pi * area_new / pow(perimeter_new, 2.0);
    out["pp_old"] = (double)-1.0 * 4.0 * pi * area_old / pow(perimeter_old, 2.0);
  }else{
    out["pp_new"] = (double)-1.0 * area_new / pow(perimeter_new, 2.0);
    out["pp_old"] = (double)-1.0 * area_old / pow(perimeter_old, 2.0);
  }
  return out;

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
		     arma::mat borderlength_mat,
		     bool discrete,
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
    aList_new = genAlConn(aList, NumericVector(new_dists.begin(), new_dists.end()));
    aList_current = genAlConn(aList, NumericVector(current_dists.begin(), current_dists.end()));
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
      psi_new += as<double>(fh_out["ssd_new"]);
      psi_old += as<double>(fh_out["ssd_old"]);

    }else if(measure == "polsby-popper"){

      List pp_out = pp_compact(new_cds, current_cds,
			       as<arma::vec>(areas_vec),
			       as<arma::vec>(boundarylist_new),
			       as<arma::vec>(boundarylist_current),
			       borderlength_mat,
			       pops,
			       aList,
			       discrete);

      // Add to psi
      psi_new += as<double>(pp_out["pp_new"]);
      psi_old += as<double>(pp_out["pp_old"]);

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
		    arma::vec county_assignments,
		    arma::vec popvec)
{

  // Initialize psi values
  double psi_new = 0.0;
  double psi_old = 0.0;

  /* We want to:
     1) Get the unique county labels
     2) For each county, loop through the unique CDs and sum over the sqrt of the share of each CD in that county
     3) Sum over counties, weight by population of that county

     4) Get the unique CD labels
     5) For each CD, loop through the unique counties and sum over the sqrt of the share of the county in that CD
     6) Sum over the CDs, weight by population of that CD
  */

  // Get unique county labels
  int i; int j; int k; int pop;
  arma::vec unique_county = unique(county_assignments);
  arma::vec pop_county(unique_county.n_elem);
  arma::uvec inds;
  for(i = 0; i < unique_county.n_elem; i++){
    pop = 0;
    inds = arma::find(county_assignments == unique_county(i));
    for(j = 0; j < inds.n_elem; j++){
      pop += popvec(inds(j));
    }
    pop_county(i) = pop;
  }

  // Loop through counties, get number of congressional districts under new and old plan
  arma::vec currentdists_incounty;
  arma::vec newdists_incounty;
  arma::vec unique_currentdists;
  arma::vec unique_newdists;
  for(i = 0; i < unique_county.n_elem; i++){
    currentdists_incounty = current_dists.elem( find(county_assignments == unique_county(i)) );
    newdists_incounty = new_dists.elem( find(county_assignments == unique_county(i)) );
    unique_currentdists = unique(currentdists_incounty);
    unique_newdists = unique(newdists_incounty);
    if(unique_currentdists.n_elem > 1){
      psi_old++;
    }
    if(unique_newdists.n_elem > 1){
      psi_new++;
    }
  }

  psi_old = (double)psi_old / unique_county.n_elem;
  psi_new = (double)psi_new / unique_county.n_elem;

  // // Get unique CD labels
  // arma::vec unique_cd = unique(new_dists);
  // arma::vec pop_cd_new(unique_cd.n_elem);
  // arma::vec pop_cd_current(unique_cd.n_elem);
  // int pop_new; int pop_current;
  // for(i = 0; i < unique_cd.n_elem; i++){
  //   pop_new = 0; pop_current = 0;
  //   inds = arma::find(new_dists == unique_cd(i));
  //   for(j = 0; j < inds.n_elem; j++){
  //     pop_new += popvec(inds(j));
  //   }
  //   inds = arma::find(current_dists == unique_cd(i));
  //   for(j = 0; j < inds.n_elem; j++){
  //     pop_current += popvec(inds(j));
  //   }
  //   pop_cd_new(i) = pop_new;
  //   pop_cd_current(i) = pop_current;
  // }

  // // Loop through the counties
  // arma::uvec county_index;
  // arma::vec pops_incounty;
  // arma::vec current_distassign_incounty;
  // arma::vec new_distassign_incounty;
  // arma::vec unique_current_dists;
  // arma::vec unique_new_dists;
  // arma::uvec inds_indistrict;

  // double ent_cd_current;
  // double ent_cd_new;
  // double ent_overcounties_current = 0.0;
  // double ent_overcounties_new = 0.0;
  // for(i = 0; i < unique_county.n_elem; i++){

  //   // Get the new and old CD assignments of the current and new plans,
  //   // and their labels
  //   county_index = arma::find(county_assignments == unique_county(i));

  //   pops_incounty = popvec.elem(county_index);

  //   current_distassign_incounty = current_dists.elem(county_index);
  //   new_distassign_incounty = new_dists.elem(county_index);

  //   unique_current_dists = unique(current_distassign_incounty);
  //   unique_new_dists = unique(new_distassign_incounty);

  //   // Get populations of the cds in the county
  //   ent_cd_current = 0.0;
  //   for(j = 0; j < unique_current_dists.n_elem; j++){
  //     inds_indistrict = arma::find(current_distassign_incounty == unique_current_dists(j));
  //     pop = 0;
  //     for(k = 0; k < inds_indistrict.n_elem; k++){
  // 	pop += pops_incounty(inds_indistrict(k));
  //     }
  //     ent_cd_current += pow((double)pop / pop_cd_current(unique_current_dists(j)), 0.5);
  //   }
  //   ent_overcounties_current += pop_county(i) * ent_cd_current;
  //   ent_cd_new = 0.0;
  //   for(j = 0; j < unique_new_dists.n_elem; j++){
  //     inds_indistrict = arma::find(new_distassign_incounty == unique_new_dists(j));
  //     pop = 0;
  //     for(k = 0; k < inds_indistrict.n_elem; k++){
  // 	pop += pops_incounty(inds_indistrict(k));
  //     }
  //     ent_cd_new += pow((double)pop / pop_cd_new(unique_new_dists(j)), 0.5);
  //   }
  //   ent_overcounties_new += pop_county(i) * ent_cd_new;

  // }

  // // Loop through the CDs
  // arma::uvec cd_index_new;
  // arma::uvec cd_index_current;
  // arma::vec pops_incd_new;
  // arma::vec pops_incd_current;
  // arma::vec current_countyassign_indist;
  // arma::vec new_countyassign_indist;
  // arma::vec unique_current_counties;
  // arma::vec unique_new_counties;

  // double ent_county_current;
  // double ent_county_new;
  // double ent_overcds_current = 0.0;
  // double ent_overcds_new = 0.0;
  // for(i = 0; i < unique_cd.n_elem; i++){

  //   // Get the indices of the the units in the new and old plans
  //   cd_index_new = arma::find(new_dists == unique_cd(i));
  //   cd_index_current = arma::find(current_dists == unique_cd(i));

  //   pops_incd_new = popvec.elem(cd_index_new);
  //   pops_incd_current = popvec.elem(cd_index_current);

  //   current_countyassign_indist = county_assignments.elem(cd_index_current);
  //   new_countyassign_indist = county_assignments.elem(cd_index_current);

  //   unique_current_counties = unique(current_countyassign_indist);
  //   unique_new_counties = unique(new_countyassign_indist);

  //   // Get populations of the counties in the cd
  //   ent_county_current = 0.0;
  //   for(j = 0; j < unique_current_counties.n_elem; j++){
  //     inds_indistrict = arma::find(current_countyassign_indist = unique_current_counties(j));
  //     pop = 0;
  //     for(k = 0; k < inds_indistrict.n_elem; k++){
  // 	pop += pops_incd_current(inds_indistrict(k));
  //     }
  //     ent_county_current += pow((double)pop / pop_county(unique_current_counties(j)), 0.5);
  //   }
  //   ent_overcds_current += pop_cd_current(i) * ent_county_current;

  //   ent_county_new = 0.0;
  //   for(j = 0; j < unique_new_counties.n_elem; j++){
  //     inds_indistrict = arma::find(new_countyassign_indist = unique_new_counties(j));
  //     pop = 0;
  //     for(k = 0; k < inds_indistrict.n_elem; k++){
  // 	pop += pops_incd_new(inds_indistrict(k));
  //     }
  //     ent_county_new += pow((double)pop / pop_county(unique_new_counties(j)), 0.5);
  //   }
  //   ent_overcds_new += pop_cd_new(i) * ent_county_new;

  // }

  // // Calculate the psis
  // psi_new = ent_overcds_new + ent_overcounties_new;
  // psi_old = ent_overcds_current + ent_overcounties_current;

  // Create return object
  List out;
  out["countysplit_new_psi"] = psi_new;
  out["countysplit_old_psi"] = psi_old;

  return out;

}

