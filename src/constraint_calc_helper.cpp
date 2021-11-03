///////////////////////////////////////////////
// Author: Ben Fifield
// Institution: Princeton University
// Date Created: 2014/12/26
// Date Last Modified: 2015/02/26
// Purpose: Contains functions to run calculate beta constraints
///////////////////////////////////////////////

// Header files
#include <RcppArmadillo.h>
#include "redist_types.h"
#include "kirchhoff.h"
#include "tree_op.h"

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

// Edges Removed Measure
List er_compact(const Graph g, arma::vec new_dists, arma::vec current_dists, int ndists){
  NumericVector er;
  int nprec = new_dists.size();
  mat districts(nprec, 2, fill::zeros);


  for(int r = 0; r < nprec; r++){
    districts(r, 0) = new_dists(r);
    districts(r, 1) = current_dists(r);
  }

  umat udistricts = conv_to<umat>::from(districts);

  er = n_removed(g, udistricts, ndists);

  List out;
  out["er_new"] = (double) er[0];
  out["er_old"] = (double) er[1];

  return out;
}

//
List log_st_compact(const Graph g, arma::vec new_dists, arma::vec current_dists,
                    arma::uvec counties, int ndists){
  NumericVector lst;
  int nprec = new_dists.size();
  mat districts(nprec, 2, fill::zeros);


  for(int r = 0; r < nprec; r++){
    districts(r, 0) = new_dists(r);
    districts(r, 1) = current_dists(r);
  }

  umat udistricts = conv_to<umat>::from(districts);

  lst = log_st_map(g, udistricts, counties + 1, ndists);

  List out;
  out["lst_new"] = (double) lst[0];
  out["lst_old"] = (double) lst[1];

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
		     // For Edges Removed & log-st
		     int ndists,
		     const Graph &g,
		     // For log-st
		     arma::vec counties,
		     // For Fryer Holden
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
  if(measure == "fryer-holden"|| measure == "polsby-popper"){
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
  } else if(measure == "edges-removed"){

    List er_out = er_compact(g, new_dists, current_dists, ndists);

    psi_new = as<double>(er_out["er_new"]);
    psi_old = as<double>(er_out["er_old"]);
  } else {

    List lst_out = log_st_compact(g, new_dists, current_dists, conv_to<uvec>::from(counties), ndists);

    psi_new = as<double>(lst_out["lst_new"]);
    psi_old = as<double>(lst_out["lst_old"]);

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

// Function to constrain by segregating a group
List calc_psivra(arma::vec current_dists,
			 arma::vec new_dists,
			 NumericVector pops,
			 NumericVector distswitch,
			 NumericVector grouppop,
			 double tgt_min,
			 double tgt_other)
{

  /* Inputs to function:
     current_dists: vector of the current cong district assignments
     new_dists: vector of the new cong district assignments
     pops: vector of district populations
     beta_vra: strength of the beta constraint
     distswitch: vector containing the old district, and the proposed new district
     grouppop: vector of subgroup district populations

  */
  // constants
  double pow_vra = 1.5;

  // Initialize psi values
  double psi_new = 0.0;
  double psi_old = 0.0;

  // Initialize denominator
  //int T = sum(pops);
  //double pAll = (double)sum(grouppop) / T;
  //double denom = (double)2 * T * pAll * (1 - pAll);

  // Loop over congressional districts
  for(int i = 0; i < distswitch.size(); i++){

    // Initialize objects
    int oldpopall = 0;
    int newpopall = 0;
    int oldpopgroup = 0;
    int newpopgroup = 0;
    arma::uvec new_cds = find(new_dists == distswitch(i));
    arma::uvec current_cds = find(current_dists == distswitch(i));

    // vra for proposed assignments
    for(int j = 0; j < new_cds.size(); j++){
      newpopall += pops(new_cds(j));
      newpopgroup += grouppop(new_cds(j));
    }

    // vra for current assignments
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
    //psi_new += (double)(newpopall * std::abs(newgroupprop - pAll));
    //psi_old += (double)(oldpopall * std::abs(oldgroupprop - pAll));

    psi_new += (double)(std::pow(std::abs(newgroupprop - tgt_min),pow_vra)*
      std::pow(std::abs(newgroupprop - tgt_other),pow_vra));
    psi_old += (double)(std::pow(std::abs(oldgroupprop - tgt_min),pow_vra)*
      std::pow(std::abs(oldgroupprop - tgt_other),pow_vra));

  }

  // Create return object
  List out;
  out["vra_new_psi"] = psi_new;
  out["vra_old_psi"] = psi_old;

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
  int i; int j; int pop;
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



// Function to constrain by segregating a group
List calc_psipartisan(arma::vec current_dists,
                 arma::vec new_dists,
                 IntegerVector rvote,
                 IntegerVector dvote,
                 std::string measure,
                 int ndists){
  double psi_new, psi_old;
  double totvote = (double)sum(rvote) + (double)sum(dvote);

  if(measure == "efficiency-gap"){
  // Step 1: Aggregate to District Level Votes
  IntegerMatrix rcounts(ndists, 2);
  IntegerMatrix dcounts(ndists, 2);

  for(int r = 0; r < current_dists.size(); r++){
    //current
    rcounts(current_dists(r), 0) += rvote(r);
    dcounts(current_dists(r), 1) += dvote(r);
    // new
    rcounts(new_dists(r), 1) += rvote(r);
    dcounts(new_dists(r), 1) += dvote(r);
  }
  // Step 2: Compute Wasted Votes
  IntegerMatrix rwaste(ndists, 2);
  IntegerMatrix dwaste(ndists, 2);
  int minwin;
  IntegerVector netwaste(2);

  for(int c = 0; c < 2; c++){
    for(int r = 0; r < ndists; r++){
      minwin = floor((dcounts(r,c) + rcounts(r,c))/2.0)+1;
      if(dcounts(r,c) > rcounts(r,c)){
        dwaste(r,c) += (dcounts(r,c) - minwin);
        rwaste(r,c) += rcounts(r,c);
      } else{
        dwaste(r,c) += dcounts(r,c);
        rwaste(r,c) += (rcounts(r,c) - minwin);
      }
    }
  }

  netwaste = colSums(dwaste) - colSums(rwaste);

  psi_old = std::abs((double)netwaste(0)/totvote);
  psi_new = std::abs((double)netwaste(1)/totvote);

} else if(measure == "proportional-representation"){
  // Step 1: Aggregate to District Level Votes
  IntegerMatrix rcounts(ndists, 2);
  IntegerMatrix dcounts(ndists, 2);

  for(int r = 0; r < current_dists.size(); r++){
    //current
    rcounts(current_dists(r)-1, 0) += rvote(r);
    dcounts(current_dists(r)-1, 1) += dvote(r);
    // new
    rcounts(new_dists(r)-1, 1) += rvote(r);
    dcounts(new_dists(r)-1, 1) += dvote(r);
  }

  // Step 2: Calculate the target Dem percent
  double target = (double)sum(dvote)/(totvote);


  // Step 3: Calculate Number of Dem Seats
  double actual_new = 0.0;
  double actual_old = 0.0;
  for(int d = 0; d < ndists; d++){
    if(dcounts(d, 0) > rcounts(d, 0)){
      actual_old += 1.0;
    }
    if(dcounts(d, 1) > rcounts(d, 1)){
      actual_new += 1.0;
    }
  }

  // Step 4: Create psi
  psi_old = actual_old/(double)ndists;
  psi_new = actual_new/(double)ndists;
}


  // Create return object
  List out;
  out["partisan_new_psi"] = std::abs(psi_new);
  out["partisan_old_psi"] = std::abs(psi_old);

  return out;
}


// Function to calculate the direct limiting constraint
List calc_psiminority(arma::vec current_dists,
                      arma::vec new_dists,
                      NumericVector pops,
                      NumericVector grouppop,
                      int ndists,
                      NumericVector minorityprop){

  int nminority = minorityprop.size();
  double psi_old = 0.0;
  double psi_new = 0.0;


  NumericVector mins_curr(ndists);
  NumericVector mins_new(ndists);
  NumericVector pops_curr(ndists);
  NumericVector pops_new(ndists);

  // Step 1: Aggregate to District Level Pops
  for(int r = 0; r < current_dists.size(); r++){
    //current
    mins_curr(current_dists(r)) += grouppop(r);
    pops_curr(current_dists(r)) += pops(r);
    // new
    mins_new(new_dists(r)) += grouppop(r);
    pops_new(new_dists(r)) += pops(r);
  }

  // Step 2: Get sort proportions
  NumericVector minprop_curr(ndists);
  NumericVector minprop_new(ndists);

  for(int i = 0; i < ndists; i++){
    minprop_curr(i) = mins_curr(i)/pops_curr(i);
    minprop_new(i) = mins_new(i)/pops_new(i);
  }

  minprop_curr.sort(TRUE);
  minprop_new.sort(TRUE);

  // Estimate psi
  for(int i = 0; i < nminority; i++){
    psi_old += std::sqrt(std::abs(minprop_curr(i) - minorityprop(i)));
    psi_new += std::sqrt(std::abs(minprop_new(i) - minorityprop(i)));
  }


  // Create return object
  List out;
  out["minority_new_psi"] = psi_new;
  out["minority_old_psi"] = psi_old;

  return out;
}

// Function to calculate the direct limiting constraint
List calc_psihinge(arma::vec current_dists,
                      arma::vec new_dists,
                      NumericVector pops,
                      NumericVector grouppop,
                      int ndists,
                      NumericVector minorityprop){

  int nminority = minorityprop.size();
  double psi_old = 0.0;
  double psi_new = 0.0;
  double tgt_old, tgt_new, diff_old, diff_new, diff_old_check, diff_new_check;


  NumericVector mins_curr(ndists);
  NumericVector mins_new(ndists);
  NumericVector pops_curr(ndists);
  NumericVector pops_new(ndists);

  // Step 1: Aggregate to District Level Pops
  for(int r = 0; r < current_dists.size(); r++){
    //current
    mins_curr(current_dists(r)) += grouppop(r);
    pops_curr(current_dists(r)) += pops(r);
    // new
    mins_new(new_dists(r)) += grouppop(r);
    pops_new(new_dists(r)) += pops(r);
  }

  // Step 2: Get sort proportions
  NumericVector minprop_curr(ndists);
  NumericVector minprop_new(ndists);

  for(int i = 0; i < ndists; i++){
    minprop_curr(i) = mins_curr(i)/pops_curr(i);
    minprop_new(i) = mins_new(i)/pops_new(i);
  }


  // Estimate psi
  for(int i = 0; i < ndists; i++){
    diff_old = 1;
    diff_new = 1;

    for(int j = 0; j < nminority; j++){
      diff_old_check = std::fabs(minorityprop(j) - minprop_curr(i));
      diff_new_check = std::fabs(minorityprop(j) - minprop_new(i));

      if(diff_new_check <= diff_new){
        diff_new = diff_new_check;
        tgt_new = minorityprop(j);
      }
      if(diff_old_check <= diff_old){
        diff_old = diff_old_check;
        tgt_old = minorityprop(j);
      }

    }

      psi_old += std::sqrt(std::max(0.0, tgt_old - minprop_curr(i)));
      psi_new += std::sqrt(std::max(0.0, tgt_new - minprop_new(i)));
  }



  // Create return object
  List out;
  out["hinge_new_psi"] = psi_new;
  out["hinge_old_psi"] = psi_old;

  return out;
}

// Function to calculate the Quadratic Population Splits constraint
List calc_psiqps(arma::vec current_dists,
              arma::vec new_dists,
              Rcpp::NumericVector pops,
              Rcpp::IntegerVector cities,
              int ndists,
              int n_city) {

  // init objects
  IntegerMatrix tally_old(ndists, n_city);
  IntegerMatrix tally_new(ndists, n_city);
  NumericMatrix pj_old(ndists, n_city - 1);
  NumericMatrix pj_new(ndists, n_city - 1);
  NumericVector j_old(n_city - 1);
  NumericVector j_new(n_city - 1);
  NumericVector sumpj_old(n_city - 1);
  NumericVector sumpj_new(n_city - 1);


  // tally populations by distr x city
  for (int i = 0; i < current_dists.size(); i ++) {
    tally_old(current_dists(i), cities(i)) += pops(i);
    tally_new(new_dists(i), cities(i)) += pops(i);
  }

  // get the district pops
  IntegerVector dist_pop_old = rowSums(tally_old);
  IntegerVector dist_pop_new = rowSums(tally_new);

  // set up proportions based on tallies
  for (int d = 0; d < ndists; d++) {
    for (int c = 1; c <  n_city; c++) {
      // distr x city proportions
      pj_old(d, c - 1) = (double) tally_old(d, c) / dist_pop_old(d);
      pj_new(d, c - 1) = (double) tally_new(d, c) / dist_pop_new(d);

      // Sum the QPS part
      sumpj_old(c - 1) += pj_old(d, c - 1) * (1.0 - pj_old(d, c - 1));
      sumpj_new(c - 1) += pj_new(d, c - 1) * (1.0 - pj_new(d, c - 1));

      // count nonzero pops
      if (tally_old(d, c) > 0) {
        j_old(c - 1) += 1;
      }
      if (tally_new(d, c) > 0) {
        j_new(c - 1) += 1;
      }
    }
  }

  // do the normalizations
  sumpj_old = sumpj_old / (double) ndists;
  sumpj_new = sumpj_new / (double) ndists;
  j_old = log(j_old);
  j_new = log(j_new);

  // calc psi
  double psi_old = sum(sumpj_old) + sum(j_old);
  double psi_new = sum(sumpj_new) + sum(j_new);

  // Create return object
  return List::create(
    _["qps_new_psi"] = psi_new,
    _["qps_old_psi"] = psi_old
  );
}
