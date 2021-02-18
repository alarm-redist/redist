///////////////////////////////////
// Author: Ben Fifield
// Institution: Princeton University
// Date Created: 2015/02/26
// Date Modified: 2015/02/26
// Purpose: Helper functions for redist package
///////////////////////////////////

// Header files
#include <RcppArmadillo.h>
#include "constraint_calc_helper.h"

using namespace Rcpp;

// Function to calculate polsby popper
// [[Rcpp::export]]
double calc_polsbypopper(arma::uvec new_cds,
			 arma::vec areas_vec,
			 arma::vec boundarylist_new,
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
  int j; int k; int l;
  if(discrete == false){
    for(j = 0; j < new_cds.n_elem; j++){
      area_new += areas_vec(new_cds(j));
    }
  }else{
    for(j = 0; j < new_cds.n_elem; j++){
      area_new += pop_vec(new_cds(j));
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

  // Get indices of the obsevations that are both in the district and
  // on the boundary
  arma::uvec boundary_inds_new = find(boundarylist_new == 1);
  arma::ivec boundary_precs_new = arma::intersect(arma::conv_to<arma::ivec>::from(new_cds), arma::conv_to<arma::ivec>::from(boundary_inds_new));

  // Loop over the indices in boundary_precs_new and check the adjacent units for
  // whether they too lie on a boundary
  arma::vec adj_precs;
  arma::uvec in_cd;
  arma::vec adj_precs_sub;

  arma::uvec adj_precs_sub_indices;
  arma::uvec find_boundarylength_in_boundarymat;
  arma::ivec loop_over_inds;
  double perimeter_new = 0.0;
  for(j = 0; j < boundary_precs_new.n_elem; j++){

    if(discrete == false){

      // Get the adjacent indices - on the boundary and
      // greater than boundary_precs[j] to avoid double counting
      adj_precs = as<arma::vec>(aList(boundary_precs_new(j)));

      // Of the adjacent precincts, which are not in the same congressional district?
      in_cd = getIn(arma::conv_to<arma::ivec>::from(adj_precs), arma::conv_to<arma::ivec>::from(new_cds));
      adj_precs_sub = adj_precs.elem( find(in_cd == false) );
      adj_precs_sub.resize(adj_precs_sub.size() + 1);
      adj_precs_sub(adj_precs_sub.size() - 1) = -1.0;

      // Which elements of boundary_mat's first column are in adj_precs_sub, and what of second column equal the actual precinct
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

  double out;
  if(discrete == false){
    out = (double)4.0 * pi * area_new / pow(perimeter_new, 2.0);
  }else{
    out = area_new / pow(perimeter_new, 2.0);
  }

  return out;

}

// Function to calculate summary statistics (package version)
// [[Rcpp::export]]
NumericVector segregationcalc(NumericMatrix distmat,
			      NumericVector grouppop,
			      NumericVector fullpop)
{

  /* Inputs to function:
     distmat: matrix of congressional districts

     grouppop: vector of the group populations

     fullpop: vector of district populations
   */

  // Vector to hold dissimilarity indices
  NumericVector diVec(distmat.ncol());

  // Population parameters
  int T = sum(fullpop);
  double P = (double)sum(grouppop) / T;

  // Calculate denominators
  double d = (double)1 / (2 * T * P * (1 - P));

  // Get the number of unique plans
  NumericVector cd1 = distmat(_,1);
  arma::vec cdVec1 = as<arma::vec> (cd1);
  arma::vec cdLabs = arma::unique(cdVec1);

  // Range to look over for cd's
  int start;
  int end = max(cd1) + 1;
  if(min(cd1) == 1){
    start = 1;
  }else{
    start = 0;
  }

  // Loop over possible plans
  for(int i = 0; i < distmat.ncol(); i++){

    // Create dissimilarity objects
    double dissim = 0;

    // Get a plan
    NumericVector cdvec = distmat(_,i);
    arma::vec cds = as<arma::vec> (cdvec);

    // Loop over congressional districts
    for(int j = start; j < end; j++){

      // Initialize counts of groups
      int tpop = 0;
      int gpop = 0;

      // Which precincts in the plan are in this cd?
      arma::uvec findCds = find(cds == j);

      // Loop over precincts
      for(int k = 0; k < findCds.size(); k++){

	// Add population counts
	tpop += fullpop(findCds(k));
	gpop += grouppop(findCds(k));

      }

      // Get district proportions
      double p = (double)gpop / tpop;

      // Add to dissimilarity index
      dissim += (double)d * tpop * std::abs(p - P);

    }

    // Add to vector
    diVec(i) = dissim;

  }

  // Return vector
  return diVec;

}

