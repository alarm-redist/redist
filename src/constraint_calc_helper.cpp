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

// Function to calculate the strength of the beta constraint for population
List calc_betapop(arma::vec current_dists,
		  arma::vec new_dists,
		  NumericVector pops,
		  double beta_population,
		  NumericVector distswitch)
{

  /* Inputs to function 
     current_dists: vector of the current cong district assignments

     new_dists: vector of the new cong district assignments

     pops: vector of district populations

     beta_population: strength of the beta constraint

     distswitch: vector containing the old district, and the proposed new district
   */

  // Calculate parity
  double parity = (double)sum(pops) / (max(current_dists) + 1);

  // Log_e(2)
  double loge2 = log(2.0);

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
    psi_new += (double)std::abs((pop_new / parity) - 1);
    psi_old += (double)std::abs((pop_old / parity) - 1);

  }

  // Calculate the ratio
  double ratio = (double)exp(beta_population * loge2 * (psi_new - psi_old));

  // Create return object
  List out;
  out["pop_ratio"] = ratio;
  out["pop_new_psi"] = psi_new;
  out["pop_old_psi"] = psi_old;

  return out;

}

// Function to calculate the strength of the beta constraint for compactness
// Fryer and Holden 2011 RPI index
List calc_betacompact(arma::vec current_dists,
		      arma::vec new_dists,
		      NumericVector pops,
		      double beta_compact,
		      NumericVector distswitch,
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

  // Log_e(2)
  double loge2 = log(2.0);

  // Loop over the congressional districts
  for(int i = 0; i < distswitch.size(); i++){

    // Initialize objects
    double ssd_new = 0.0;
    double ssd_old = 0.0;
    arma::uvec new_cds = find(new_dists == distswitch(i));
    arma::uvec current_cds = find(current_dists == distswitch(i));

    // SSD for new partition
    for(int j = 0; j < new_cds.size(); j++){
      for(int k = j + 1; k < new_cds.size(); k++){
	ssd_new += (double)ssdmat(new_cds(j),new_cds(k)) *
	  pops(new_cds(j)) * pops(new_cds(k));
      }
    }

    // SSD for old partition
    for(int j = 0; j < current_cds.size(); j++){
      for(int k = j + 1; k < current_cds.size(); k++){
	ssd_old += (double)ssdmat(current_cds(j),current_cds(k)) *
	  pops(current_cds(j)) * pops(current_cds(k));
      }
    }

    // Add to psi
    psi_new += ssd_new;
    psi_old += ssd_old;

  }

  // Normalize psi
  psi_new = (double)psi_new / denominator;
  psi_old = (double)psi_new / denominator;

  // Calculate ratio
  double ratio = (double)exp(beta_compact * loge2 * (psi_new - psi_old));

  // Create return object
  List out;
  out["compact_ratio"] = ratio;
  out["compact_new_psi"] = psi_new;
  out["compact_old_psi"] = psi_old;

  return out;

}

// Function to constrain by segregating a group
List calc_betasegregation(arma::vec current_dists,
			  arma::vec new_dists,
			  NumericVector pops,
			  double beta_segregation,
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

  // Log_e(2)
  double loge2 = log(2.0);

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

  // Get mh ratio
  double ratio = (double)exp(beta_segregation * loge2 * (psi_new - psi_old));

  // Create return object
  List out;
  out["segregation_ratio"] = ratio;
  out["segregation_new_psi"] = psi_new;
  out["segregation_old_psi"] = psi_old;

  return out;

}

// Function to constrain on plan similarity to original plan
List calc_betasimilar(arma::vec current_dists,
		      arma::vec new_dists,
		      arma::vec orig_dists,
		      double beta_similar,
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

  // Log_e(2)
  double loge2 = log(2.0);

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

  // Get MH ratio
  double ratio = (double)exp(beta_similar * loge2 * (psi_new - psi_old));

  // Rcout << "New psi" << psi_new << std::endl;
  // Rcout << "Old psi" << psi_old << std::endl;
  // Create return object
  List out;
  out["similar_ratio"] = ratio;
  out["similar_new_psi"] = psi_new;
  out["similar_old_psi"] = psi_old;

  return out;
}

