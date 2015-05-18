/////////////////////////////////////
// Author: Ben Fifield
// Institution: Princeton University
// Date Created: 2014/12/17
// Date Modified: 2015/02/26
// Purpose: Supporting functions for swMH() function in redist
/////////////////////////////////////

// Header files
#include <RcppArmadillo.h>
#include "make_swaps_helper.h"
#include "constraint_calc_helper.h"

using namespace Rcpp;

// Function to generate initial vector of populations
NumericVector init_pop(NumericVector popvec,
		       arma::vec cds)
{

  /* Inputs to function:
     cds: Vector of congressional district populations

     popvec: Vector of populations
   */ 

  // Get number of cds
  int ncds = cds.max() + 1;

  // Create container vector
  NumericVector distpop(ncds);

  // Loop through cd assignments
  for(int i = 0; i < ncds; i++){

    // Initialize population count
    int pop = 0;

    // Get indices of cds 
    arma::uvec cd_i_ind = find(cds == i);
    
    // Loop through cd_i_ind, get population values
    for(int j = 0; j < cd_i_ind.n_elem; j++){
      pop += popvec(cd_i_ind(j));
    }

    // Put in distpop
    distpop(i) = pop;

  }

  return distpop;

}

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

  // Loop through precincts
  for(int i = 0; i < cds.size(); i++){

    // For precinct i, get adjacent precincts
    NumericVector avec = aList(i);
    
    // Get precinct i's congressional district
    int cd_i = cds(i);

    // Initialize empty vector
    NumericVector avec_cd;

    // Loop through avec to identify which are in same cd
    for(int j = 0; j < avec.size(); j++){
      
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
NumericVector findBoundary(List fullList,
			   List conList)
{

  /* Inputs to function:
     fullList: Full adjacency list of geographic units

     conList: Adjacency list of geographic units within cong district
   */

  // Initialize container vector of 0's (not boundary) and 1's (boundary)
  NumericVector isBoundary(fullList.size());

  // Loop through aList
  for(int i = 0; i < fullList.size(); i++){

    // Get vectors of full and cd-connected components for precinct i
    NumericVector full = fullList(i);
    NumericVector conn = conList(i);

    // Compare lengths - if conn < full, then boundary unit
    if(full.size() > conn.size()){
      isBoundary(i) = 1;
    }
    
  }

  return isBoundary;

}

// Function to make unidirectional adjacency list bidirectional
List add_ties(List aList){
  
  // Loop through vectors in aList
  for(int i = 0; i < aList.size(); i++){

    // Get i'th entry in list
    NumericVector list1 = aList(i);
    
    // Loop through elements in list1
    for(int j = 0; j < list1.size(); j++){

      // Extract adjacency vector for j'th element of i's adjacency list
      NumericVector list2 = aList(list1(j));
      
      // Check if list 2 includes i
      if(is_true(any(list2 == i)) == FALSE){

	// If not included, add to adjacency vector
	list2.push_back(i);

	// Modify aList to include new adjacency vector
	aList(list1(j)) = list2;

      }
    
    }
 
  }

  return aList;

}

// Function to cut edges of adjacency list probabilistically - Step 2 of swMH
List cut_edges(List aList_con,
	       double eprob)
{

  /* Inputs to function:
     aList_con: adjacency list within cong district

     eprob: edgecut probability (transformed into 1-eprob in function)
   */

  // Create threshold
  double threshold_prob = 1 - eprob;

  // Define lists to store cut-edge and uncut-edge vectors
  List aList_uncut(aList_con.size());
  List aList_cut(aList_con.size());

  // Define list to store output of both lists

  // Loop through elements of aList_con
  for(int i = 0; i < aList_con.size(); i++){

    // Extract i'th vector in list
    NumericVector cc_vec_i_all = aList_con(i);

    // Subset cc_vec_i to elements > i
    NumericVector cc_vec_i = cc_vec_i_all[cc_vec_i_all > i];

    // For each element in vector, take random draw from [0,1] uniform
    arma::vec draws = runif(cc_vec_i.size());

    // Create container vectors of cut and uncut edges
    NumericVector cut;
    NumericVector uncut;

    // Loop through elements of cc_vec_i and compare to entry in draws
    for(int j = 0; j < cc_vec_i.size(); j++){
      
      // Compare to threshold_prob - if draws < thresh, cut edge, else uncut
      if(draws(j) < threshold_prob){
	cut.push_back(cc_vec_i(j));
      } else{
	uncut.push_back(cc_vec_i(j));
      }

    }

    /* Here - look at lines 1201-1212 in original code. Modifying original
       alConnected to remove edges that are cut, but isn't this just the 
       uncut list (which will be aList_postcut? Skipping this bit for now */
    
    // Store vectors in container lists
    aList_uncut(i) = uncut;
    aList_cut(i) = cut;

  }
  
  // Add ties to aList_uncut, aList_cut
  List aList_uncut_bd = add_ties(aList_uncut);
  List aList_cut_bd = add_ties(aList_cut);

  // Return contents
  List out;
  out["connectedlist"] = aList_uncut_bd;
  out["cutedgelist"] = aList_cut_bd;
  
  return out;

}

/* Function to run breadth-first search, returning only sets of connected 
   components that reside on the boundary of the districts */
List bsearch_boundary(List aList,
		      arma::vec boundary)
{

  /* Inputs to function:
     aList: adjacency list

     boundary: vector of boundary element indicators (as arma)
   */

  // Get indices of boundary units
  arma::uvec boundary_indices = find(boundary == 1);

  // Container - outputted of breadth search, a list
  List bsearch;

  // Container - partition vector, gets added to bsearch when queue is empty
  NumericVector partition;

  // Set mark vector - ledger of which indices have been reached
  NumericVector mark(aList.size());

  // Set queue vector
  NumericVector q;

  // Initialize breadth search with first element in boundary_indices
  mark(boundary_indices(0)) = boundary_indices(0);
  partition.push_back(boundary_indices(0));
  q = aList(boundary_indices(0));

  // Begin do{} loop - run until number of elements in boundary_indices is 0
  do{

    // Begin while{} loop - run until q is empty
    while(q.size() > 0){
      
      // Dequeue first element in queue
      int u = q(0);

      // Mark that element in ledger
      mark(u) = u;

      // Check if element is in the partition - add to partition if false
      bool in_part = is_true(any(partition == u));
      if(in_part == FALSE){
	partition.push_back(u);
      }
      
      // Get adjacency vector for unit u
      NumericVector adj_u = aList(u);

      // Loop through elements of adj_u, add to queue and mark if not reached
      if(adj_u.size() > 0){
	
	// Start loop
	for(int i = 0; i < adj_u.size(); i++){
	  
	  // Reach element v
	  int v = adj_u(i);

	  /* Check if already reached - if false, mark, add to partition, and
	     add to queue */
	  if(is_true(any(mark == v)) == FALSE){
	    mark(v) = v;
	    partition.push_back(v);
	    q.push_back(v);
	  }

	}

      }

      // Erase dequeued element from queue when done searching
      q.erase(q.begin());

    }

    // Handling an empty queue
    if(q.size() == 0){

      /* First, find boundary units that are in the reached partition and
	 remove them from boundary_units vector */
      for(int i = boundary_indices.n_elem - 1; i >= 0; i--){
	if(is_true(any(partition == boundary_indices(i))) == TRUE){
	  boundary_indices.shed_row(i);
	}
      }
      
      // Store the partition, clear partition vector
      bsearch.push_back(partition);
      partition.erase(partition.begin(), partition.end());

      // Re-initialize breadth search from new starting value if nonempty
      if(boundary_indices.n_elem > 0){
	q = aList(boundary_indices(0));
	mark(boundary_indices(0)) = boundary_indices(0);
	partition.push_back(boundary_indices(0));
      }

    }

  }while(boundary_indices.n_elem > 0);

  List out;
  out["bsearch"] = bsearch;
  out["npartitions"] = bsearch.size();

  return out;

}

/* Function to draw p for the number of connected components */
int draw_p(int lambda)
{

  /* Inputs to function:
     lambda: lambda parameter
   */

  int p;
  if(lambda > 0){
    p = R::rpois(lambda);
    p++;
  } else{
    p = 1;
  }

  return p;

}

/* Function to draw p separate, noncontiguous connected components from the
   output of the boundary breadth search. These are candidate swaps
   to form the next iteration of the markov chain. Function returns
   the proposed district assignments that will be accepted or rejected */
List make_swaps(List boundary_cc,
		List aList, 
		NumericVector cds_old,
		NumericVector cds_orig,
		NumericVector pop_vec, 
		NumericVector cd_pop_vec,
		NumericVector group_pop_vec,
		NumericMatrix ssdmat,
		double minparity,
		double maxparity, 
		int p, 
		double eprob,
		double beta_population,
		double beta_compact,
		double beta_segregation,
		double beta_similar,
		double ssd_denominator)
{

  /* Inputs to function:
     boundary_cc: Connected components on district boundaries

     aList: full adjacency list

     cds_old: Current cong district assignments

     cds_orig: original cong district assignments. For similarity constraint

     pop_vec: unit populations

     cd_pop_vec: congressional district populations

     group_pop_vec: populations of groups in geographic units

     ssdmat: sum of squared distance matrix

     minparity, maxparity: population parity - min and max

     p: parameter p for the number of swaps

     eprob: edgecut probability

     beta_population: strength of constraint for achieving population parity.

     beta_compact: strength of constraint for achieving compactness

     beta_segregation: strength of constraint for segregating subgroup

     beta_similar: strength of constraint for similarity to orig plan

     ssd_denominator: normalizing constant for sum of squared distance psi

   */
  
  // Initialize objects for swap //
  NumericVector cds_prop = clone(cds_old);
  NumericVector cdspop_prop = clone(cd_pop_vec);
  NumericVector accepted_partitions;
  NumericVector prop_cd_pops;
  NumericVector cds_test;

  // Initialize metropolis-hastings probabilities
  double mh_prob = 1.0;

  // Number of unique congressional districts
  int ndists = max(cds_old) + 1;

  // Break indicators
  int breakp = 0; 
  int numaccept = 0;
  int goodprop = 0;
  
  // Begin loop over p
  for(int i = 0; i < p; i++){

    // While there are still possible components available
    int breakwhile = 0;
    int curr_cd;
    int prop_cd;
    NumericVector prop_partitions;
    while(boundary_cc.size() > 0){
      
      // (1) - select a connected component from boundary_cc randomly
      arma::vec rand_sample_index = runif(1, 0, 1000000000);
      int sample_index = fmod(rand_sample_index(0), boundary_cc.size());

      prop_partitions = boundary_cc(sample_index);
      boundary_cc.erase(sample_index);
      curr_cd = cds_prop(prop_partitions(0));
      
      /* (2) - check to see if that connected component is adjacent to one
	 already selected. Also, gathers adjacent congressional districts */
      List adjcheck_out = adjcheck_propcd(aList, 
					  prop_partitions,
					  accepted_partitions,
					  cds_prop);
      
      // If adjacent to already-accepted unit, then invalid - return to top
      int adjcheck = as<int>(adjcheck_out["adjacency_check"]);
      if(adjcheck == 1){
	continue;
      }

      /* (3) - Check to see if the proposed swap eliminates its old
	 congressional district, conditional on valid previous swaps. 
	 If invalid, return to the top */
      int elimcheck = elim_check(prop_partitions, cds_prop);
      if(elimcheck == 1){
	continue;
      }

      /* (4) - check to see if the proposed swap splits its old 
	 congressional district - make arbitrary cd assignment just to 
	 test (this is not the sample stage) */
      // Preprocess by arbitrarily assigning prop_partition to new cd
      NumericVector possible_cd_swaps = as<NumericVector>(adjcheck_out["proposed_cds"]);
      NumericVector cds_splittest = clone(cds_prop);
      for(int j = 0; j < prop_partitions.size(); j++){
	cds_splittest(prop_partitions(j)) = possible_cd_swaps(0);
      }
      // Get adjacency list
      List aList_testsplit = genAlConn(aList, cds_splittest);
      // Get number of connected components
      int num_cds = countpartitions(aList_testsplit);
      if(num_cds != ndists){
	continue;
      }
      
      // (5) - propose to swap into a new (adjacent) congressional district
      int numcds_test = possible_cd_swaps.size();
      
      // Loop over elements in propcds - try each one
      for(int j = 0; j < numcds_test; j++){

	// Draw an element from possible_cds_swaps
	if(possible_cd_swaps.size() > 1){
	  arma::vec rand_test_cd_ind = runif(1, 0, 1000000000);
	  int test_cd_ind = fmod(rand_test_cd_ind(0), possible_cd_swaps.size());
	  
	  prop_cd = possible_cd_swaps(test_cd_ind);
	  possible_cd_swaps.erase(test_cd_ind);
	} else{
	  prop_cd = possible_cd_swaps(0);
	  possible_cd_swaps.erase(0);
	}
	
	// Create a test cd vector - change cds of test partition to prop_cd
	cds_test = clone(cds_prop);
	for(int k = 0; k < prop_partitions.size(); k++){
	  cds_test(prop_partitions(k)) = prop_cd;
	}

	/* (6) - Check to see if proposed swap would violate the imposed
	   population constraint */
	prop_cd_pops = update_distpop(prop_partitions, 
				      pop_vec, 
				      prop_cd, 
				      curr_cd,
				      cdspop_prop);

	int paritycheck = 0;
	int k;
	for(k = 0; k < prop_cd_pops.size(); k++){
	  if((prop_cd_pops(k) >= minparity) && (prop_cd_pops(k) <= maxparity)){
	    paritycheck++;
	  } 
	}
	if(paritycheck == prop_cd_pops.size()){
	  breakwhile++;
	  numaccept++;
	  break;
	}
	
      }

      // Check for if we have enough valid swaps
      if(numaccept == p){
	goodprop++;
	breakp++;
	break;
      }
      
      // If good proposal, but need more
      if(breakwhile == 1){
	break;
      }
      
    }

    // If we run out of good partitions
    if((boundary_cc.size() == 0) & (numaccept < p)){
      break;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////
    // Accept if (2-6) are satisfied. First, update mh_prob and all constraints. Then, //
    // change cds_prop, cd_pop_vec, and add to accepted_partitions //////////////////////
    /////////////////////////////////////////////////////////////////////////////////////
    
    // Create vector of the cong district swap. First entry is old district, second entry
    // is proposed district
    NumericVector cd_pair(2);
    cd_pair(0) = curr_cd;
    cd_pair(1) = prop_cd;
    
    // Update metropolis-hastings probabilities
    mh_prob = update_mhprob(prop_partitions,
			    aList,
			    cds_prop,
			    prop_cd,
			    eprob,
			    mh_prob);
    
    // Calculate beta constraints
    double population_constraint = 1.0;
    if(beta_population != 0.0){
      population_constraint = as<double>(calc_betapop(cds_prop,
						      cds_test,
						      pop_vec,
						      beta_population,
						      cd_pair)["pop_ratio"]);
    }
    double compact_constraint = 1.0;
    if(beta_compact != 0.0){
      compact_constraint = as<double>(calc_betacompact(cds_prop,
						       cds_test,
						       pop_vec,
						       beta_compact,
						       cd_pair,
						       ssdmat,
						       ssd_denominator)["compact_ratio"]);
    }
    double segregation_constraint = 1.0;
    if(beta_segregation != 0.0){
      segregation_constraint = as<double>(calc_betasegregation(cds_prop,
							       cds_test,
							       pop_vec,
							       beta_segregation,
							       cd_pair,
							       group_pop_vec)["segregation_ratio"]);
    }
    double similar_constraint = 1.0;
    if(beta_similar != 0.0){
      similar_constraint = as<double>(calc_betasimilar(cds_prop,
						       cds_test,
						       cds_orig,
						       beta_similar,
						       cd_pair)["similar_ratio"]);
    }
    
    // Multiply mh_prob by constraint values
    mh_prob = (double)mh_prob * population_constraint * compact_constraint *
      segregation_constraint * similar_constraint;
    
    // Update cd assignments and cd populations
    cds_prop = cds_test;
    cdspop_prop = prop_cd_pops;
    
    // Push back prop_partition to accepted_partitions
    for(int j = 0; j < prop_partitions.size(); j++){
      accepted_partitions.push_back(prop_partitions(j));
    }

    // If we get enough partitions
    if(breakp == 1){
      break;
    }
    
  }
  
  // Create returned list
  List out;
  out["proposed_partition"] = cds_prop;
  out["mh_prob"] = mh_prob;
  out["updated_cd_pops"] = cdspop_prop;
  out["goodprop"] = goodprop;
  
  return out;
  
}

// Function to accept or reject swaps
int mh_decision(double mh_prob)
{

  /* Inputs to function:
     mh_prob: metropolis-hastings probability
  */
  
  // Initialize decision
  int decision = 0;

  // Get acceptance probability
  double acc_prob;
  if(mh_prob < 1){
    acc_prob = mh_prob;
  } else{
    acc_prob = 1;
  }

  // Draw from uniform
  arma::vec draw_prob = runif(1);

  // Make decision
  if(draw_prob(0) <= acc_prob){
    decision++;
  }

  return decision;

}

// Function that applies the Geyer Thompson algorithm for simulated tempering
List changeBeta(arma::vec betavec,
		double beta,
		double constraint,
		NumericVector weights,
		int adjswap = 1)
{
  
  /* Inputs to function 
     betavec: vector of possible betas

     beta: current value of the beta constraint

     constraint: the evaluation of the constraint on the current plan

     weights: priors on the betas

     adjswap: flag - do we want adjacent swaps? default to 1
   */
  
  // Find beta in betavec
  arma::uvec findBetaVec = find(betavec == beta);
  int findBeta = findBetaVec(0);

  // Object to test whether beta is at RHS of vector
  int betaLoc = betavec.size() - 1;

  // Get transition probabilities and propose a new beta
  double qij;
  double qji;
  double wi;
  double wj;
  double propBeta;

  // Procedure if conducting adjacent swaps
  if(adjswap == 1){
    if(findBeta == 0){ // At first element in betavec
      qij = 1;
      qji = .5;
      wi = weights(0);
      wj = weights(1);
      propBeta = betavec(1);
    } else if(findBeta == betaLoc){ // At last element in betavec
      qij = 1;
      qji = .5;
      wi = weights(betaLoc);
      wj = weights(betaLoc - 1);
      propBeta = betavec(betaLoc - 1);
    } else{ // Anywhere in the middle of betavec
      qij = .5;
      qji = .5;
      wi = weights(findBeta);
      arma::vec betaswitch = runif(1);
      if(betaswitch(0) < .5){
	propBeta = betavec(findBeta - 1);
	wj = weights(findBeta - 1);
      }
      if(betaswitch(0) >= .5){
	propBeta = betavec(findBeta + 1);
	wj = weights(findBeta + 1);
      }
    }
  } else{
    // Procedure if not conducting adjacent swaps
    // qij = qji in non-adjacent framework, don't have to worry abt end units
    qij = 1;
    qji = 1;

    // Draw element from betavec
    arma::vec rand_randindex = runif(1, 0, 1000000000);
    int randindex = fmod(rand_randindex(0), betaLoc);

    // Weight wi 
    wi = weights(findBeta);

    // Draw the proposed beta value
    if(randindex < findBeta){
      propBeta = betavec(randindex);
      wj = weights(randindex);
    } else{
      propBeta = betavec(randindex + 1);
      wj = weights(randindex + 1);
    }

  }

  // Accept or reject the proposal
  double mhprobGT = (double)exp(constraint * (propBeta - beta)) * wj / wi * qji / qij;
  if(mhprobGT > 1){
    mhprobGT = 1;
  }
  arma::vec testkeepGT = runif(1);
  int decision = 0;
  if(testkeepGT(0) <= mhprobGT){
    decision++;
    beta = propBeta;
  }

  // Create output
  List out;
  out["beta"] = beta;
  out["mh_decision"] = decision;
  out["mh_prob"] = mhprobGT;

  return out;

}

