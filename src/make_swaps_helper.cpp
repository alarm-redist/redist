///////////////////////////////////////////////
// Author: Ben Fifield
// Institution: Princeton University
// Date Created: 2014/12/26
// Date Last Modified: 2015/02/26
// Purpose: Contains functions to run make_swaps function in sw_mh_helper.cpp
/////////////////////////////////////////////// 

// Header files
#include <RcppArmadillo.h>
#include "constraint_calc_helper.h"

using namespace Rcpp;

/* Function to check adjacency of a randomly selected connected component
   against connected components already selected. Also gets candidate
   congressional district swaps if valid */
List adjcheck_propcd(List aList,
		     NumericVector prop_partitions,
		     NumericVector accepted_partitions,
		     NumericVector cds)
{  
  /* Inputs to function:
     aList: Full adjacency list

     prop_partitions: The proposed partition to be swapped

     accepted_partitions: Vector of district ID's that have been acccepted

     cds: Vector of cd assignments - has to be cds_prop to 
     avoid complications with not recognizing splitting, elimination
   */
  
  // Initialize adjacency check value
  int adj_check = 0;

  // Initialize vector of proposed congressional cds
  NumericVector prop_cds;
  
  // Get current cd of the proposed partition
  int current_cd = cds(prop_partitions(0));

  /* Loop over units in prop_partitions, test to see if any are adjacent to 
     any units in accepted_partitions */
  for(int i = 0; i < prop_partitions.size(); i++){

    // For i'th unit in prop_partitions, get the adjacency list
    NumericVector adj_units = aList(prop_partitions(i));

    // Loop to see if any of the indices in adj_units are in accepted_partitions
    for(int j = 0; j < adj_units.size(); j++){

      // See if element j of adj_units is in accepted_partitions
      bool test_adj = is_true(any(accepted_partitions == adj_units(j)));

      /* If true, iterate adj_check to 1 and break the loop to throw out 
	 the partition */
      if(test_adj == TRUE){
	adj_check++;
	break;
      }

      // If not true, then look at the congressional districts of adjacent units
      // Is unit j's congressional district equal to current_cd?
      bool same_cd = cds(adj_units(j)) == current_cd;
      // Would this be a new addition to prop_cds?
      bool new_cd = is_true(any(prop_cds == cds(adj_units(j))));

      // If both conditions are false, add to prop_cds
      if((same_cd == FALSE) && (new_cd == FALSE)){
	prop_cds.push_back(cds(adj_units(j)));
      }

    }

    // If partition is found to be adjacent, do not test any more
    if(adj_check == 1){
      break;
    }

  }

  // Create output from function
  List out;
  out["adjacency_check"] = adj_check;
  out["proposed_cds"] = prop_cds;

  return out;

}

// Function to do the elimination check
int elim_check(NumericVector prop_partition,
	       NumericVector cds)
{

  /* Inputs to function:
     prop_partition: Proposed partition

     cds: Vector of congressional districts - 
     using accepted partitions
   */ 
  
  // Indicator for elimimation
  int elimcheck = 0;

  // Get current congressional district
  int current_cd = cds(prop_partition(0));

  // Get length of cds that is of that cd
  NumericVector subcd = cds[cds == current_cd];
  int current_cd_size = subcd.size();

  // If prop_partition size is equal to number of units of that cd, then elim.
  if(current_cd_size == prop_partition.size()){
    elimcheck++;
  }

  return elimcheck;

}

// Function to generate adjacency graph and count clusters
// [[Rcpp::export]]
int countpartitions(List aList) 
{   

  //Takes an adjacency list,
  //The vector of subset nodes
  //The number of subset nodes
						
  //initialize connCompVec   
  //Initialize visited indices
  IntegerVector visitedInd(aList.size());
  int indexVisit = 0;
  
  //Initialize connected components
  IntegerVector currConnComp(aList.size());

  //Initialize the number of connected components
  int numConnComp = 0;
  
  //Loop over nodes
  for(int i = 0; i < aList.size(); i++){
    
    //If i has not been visited...
    if(visitedInd[i] == 0){
      
      //List i as visited
      visitedInd[i] = 1;

      //Increase the number of connected components
      numConnComp++;

      //Add i to the connected component list
      currConnComp[indexVisit] = i;
      
      //increase index visit
      indexVisit++;
      
      //Count the number of nodes in the current connected component
      int nodeCount = indexVisit - 1;
      
      //Initialize a stopping variable:
      int toStop = 0;

      //While we don't stop
      while(toStop == 0){
	
	//get the neighbors of the next current comp
	IntegerVector listNeighs = aList[currConnComp[nodeCount]];
	
	//If listNeighs does not have length zero...
	int listLength = listNeighs.size();
	if(listLength > 0){
	  
	  //Add nodes of listLength to currConnComp
	  //and mark nodes as visited
	  for(int j = 0; j < listLength; j++){
	    if( visitedInd[listNeighs[j]] == 0){
	      currConnComp[indexVisit] = listNeighs[j];
	      visitedInd[listNeighs[j]] = 1;

	      //Increment indexVisit
	      indexVisit++;
	    }
	  }
	}
	
	//Increment nodeCount
	nodeCount++;

	//If currConnComp[nodeCount] is zero, then we must have new connected component
	//Also stop if we have too many guys.
	if(nodeCount == aList.size()){
	  toStop = 1;
	}
	else if(currConnComp[nodeCount] == 0 ){
	  toStop = 1;
	}
      }
    }
  }
  
  return numConnComp;
  
}

// Function to update district populations
NumericVector update_distpop(NumericVector prop_partition,
			     NumericVector unitpop_vec, 
			     int prop_cd,
			     int curr_cd,
			     NumericVector distpop_vec)
{

  /* Inputs to function:
     prop_partition: Proposed partition to be swapped

     unitpop_vec: population vector for units

     prop_cd: proposed cong district for prop_partition

     curr_cd: old cong district for prop_partition

     distpop_vec: Vector of cong district populations
   */

  // Clone distpop_vec
  NumericVector distpop_vec_clone = clone(distpop_vec);
  
  // Current population, proposed district population
  int currpop = distpop_vec_clone(curr_cd);
  int proppop = distpop_vec_clone(prop_cd);

  // Loop through prop_partition
  for(int i = 0; i < prop_partition.size(); i++){
    currpop -= unitpop_vec(prop_partition(i));
    proppop += unitpop_vec(prop_partition(i));
  }

  // Put back in distpop_vec
  distpop_vec_clone(curr_cd) = currpop;
  distpop_vec_clone(prop_cd) = proppop;

  return distpop_vec_clone;

}

// Function to update the metropolis-hastings probability for a swap
double update_mhprob(NumericVector prop_partition,
		     List aList,
		     arma::vec cds,		     
		     int prop_cd,
		     double eprob, 
		     double mh_prob)
{

  /* Inputs to function:
     prop_partition: Proposed partition to swap

     aList: Full adjacency list

     cds: original vec of cong districts

     prop_cd: proposed cong district

     eprob: edgecut probability

     mh_prob: current metropolis-hastings probability
   */

  // Initialize c1, c2
  int c1 = 0;
  int c2 = 0;
	
  // Loop through prop_partition
  for(int i = 0; i < prop_partition.size(); i++){

    // Get adjacency vector
    NumericVector adj_vec = aList(prop_partition(i));

    // Loop throgh elements of adj_vec
    for(int j = 0; j < adj_vec.size(); j++){

      /* Calculate C(V_0, V_l' \ V_0) - add 1 if adjacent to switched
	 partition, and your old cd assignment is same as proposed switch */
      if(cds(adj_vec(j)) == prop_cd){
	c1++;
      }

      /* Calculate C(V_0, V_l \ V_0) - add 1 if you are adjacent to the 
	 proposed switch, if your cd assignment is the same as the old cong
	 district, and you are not in the switch partition */
      if((cds(adj_vec(j)) == cds(prop_partition(0))) &&
	 (is_true(any(prop_partition == adj_vec(j))) == FALSE)){
	c2++;
      }

    }

  }

  // Recalculate mh probability
  mh_prob = (double)mh_prob * ((double)pow(1 - eprob, c1) / pow(1 - eprob, c2));

  return mh_prob;

}

// Function to calculate squared distance matrix (written by Jonathan Olmsted, NPD Group)
// [[Rcpp::export]]
NumericMatrix calcPWDh (NumericMatrix x)
{

  int nrows = x.nrow() ;
  int ncols = x.nrow() ;
  NumericMatrix out(nrows, ncols) ;
  double rad = 3963.1676 ;
  double pi = 3.141592653589793238463 ;
  for(int arow = 0; arow < nrows; arow++) {
    for(int acol = 0; acol < ncols; acol++) {
      double phi1 = x(arow, 0) * pi / 180 ;
      double phi2 = x(acol, 0) * pi / 180 ;
      double lambda1 = x(arow, 1) * pi / 180 ;
      double lambda2 = x(acol, 1) * pi / 180;
      double q1 = 2 * rad ;
      double q2 = pow(sin((phi1 - phi2) / 2), 2) ;
      double q3 = pow(sin((lambda1 - lambda2) / 2), 2) ;
      double q4 = cos(phi1) * cos(phi2) ;
      out(arow, acol) = q1 * asin(sqrt(q2 + q4 * q3)) ;
    }
  }
  
  return out;

}

// Function to calculate deviation from original cd vector
NumericVector diff_origcds(NumericMatrix mat,
			   NumericVector cds){

  /* Inputs to function:
     mat: Matrix of congressional district assignments

     cds: Vector of congressional district assignments to serve as baseline
   */

  // Get length of cd vector
  unsigned int len_cds = cds.size();

  // Convert cds to arma
  arma::uvec cds_arma = as<arma::uvec>(cds);

  // Initialize objects
  unsigned int i; unsigned int j; unsigned int k = mat.ncol(); arma::vec plan;
  arma::uvec compare; NumericVector store_compare(k);
  
  // Start loop over mat columns
  for(i = 0; i < k; i++){

    // Get the plan
    plan = mat(_,i);

    // Compare matrices
    compare = (plan == cds_arma);

    // Sum up and divide by len_cds
    store_compare(i) = (double)sum(compare) / len_cds;
      
  }

  return store_compare;

}

// Function to calculate deviation from parity from cd matrix
NumericVector distParity(NumericMatrix mat,
			 NumericVector popvec)
{

  /* Inputs to function:
     mat: Matrix of congressional district assignments

     popvec: vector of geographic unit populations

   */

  // Get the unique cd labels
  NumericVector labs = unique(mat(_,0));

  // Calculate parity
  double parity = (double)sum(popvec) / labs.size();

  // Container vector of plan deviations
  NumericVector plandevs(mat.ncol());

  // Loop through plans
  for(int i = 0; i < mat.ncol(); i++){

    // Get plan
    arma::vec plan = mat(_,i);

    // Loop through assignments
    double maxdev = 0.0;
    for(int j = 0; j < labs.size(); j++){

      arma::uvec assignments = find(plan == labs(j));

      // Loop over precincts in plan
      int distpop = 0;
      for(int k = 0; k < assignments.size(); k++){

	distpop += popvec(assignments(k));

      }

      // Get deviation from parity
      double plandev = std::abs((double)((double)distpop - parity) / parity);
      if(plandev > maxdev){
	maxdev = plandev;
      }

    }

    // Store maxdev in plandevs
    plandevs(i) = maxdev;

  }

  return plandevs;

}

