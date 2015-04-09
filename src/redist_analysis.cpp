///////////////////////////////////
// Author: Ben Fifield
// Institution: Princeton University
// Date Created: 2015/02/26
// Date Modified: 2015/02/26
// Purpose: Helper functions for redist package
///////////////////////////////////

// Header files
#include <RcppArmadillo.h>

using namespace Rcpp;

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

