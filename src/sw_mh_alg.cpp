/////////////////////////////////////
// Author: Ben Fifield
// Institution: Princeton University
// Date Created: 2014/12/15
// Date Modified: 2015/02/26
// Purpose: Code for swMH() function in redist package
/////////////////////////////////////

// Header files
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <time.h>
#include <R.h>
#include "sw_mh_helper.h"
#include "make_swaps_helper.h"
#include "constraint_calc_helper.h"
#include "redist_analysis.h"

using namespace Rcpp;

/* Primary function to run redistricting algorithm. An implementation of
   Algorithm 1 in Barbu and Zhu (2005) */
// [[Rcpp::export]]
List swMH(List aList,
	  NumericVector cdvec,
	  NumericVector cdorigvec,
	  NumericVector popvec,
	  NumericVector grouppopvec,
	  int nsims,
	  double eprob,
	  double pct_dist_parity,
	  NumericVector beta_sequence,
	  NumericVector beta_weights,
	  NumericMatrix ssdmat,
	  int lambda = 0,
	  double beta_population = 0.0,
	  double beta_compact = 0.0,
	  double beta_segregation = 0.0,
	  double beta_similar = 0.0,
	  int anneal_beta_population = 0,
	  int anneal_beta_compact = 0,
	  int anneal_beta_segregation = 0,
	  int anneal_beta_similar = 0,
		std::string adapt_beta = "none",
	  int adjswap = 1,
	  int exact_mh = 0,
	  int adapt_eprob = 0,
	  int adapt_lambda = 0,
		int num_hot_steps = 0,
	  int num_annealing_steps = 0,
	  int num_cold_steps = 0)
{

  /* Inputs to function:
     aList: adjacency list of geographic units

     cdvec: initial vector of congressional district assignments

     popvec: vector of populations for each geographic unit

     grouppopvec: vector of subgroup populations for each geographic unit

     nsims: number of simulations

     eprob: edgecut probability

     pct_dist_parity: strength of population parity requirement

     beta_sequence: sequence of betas to anneal ove

     beta_weights: prior weights on the beta sequence

     ssdmat: matrix of squared distances between geographic units. For constraining
     on compactness

     lambda: parameter for conducting multiple swaps

     beta_population: strength of constraint for achieving population parity.

     beta_compact: strength of constraint for achieving district compactness

     beta_segregation: strength of constraint for packing group into district

     beta_similar: strength of constraint for examining plans similar to original district

     anneal_beta_population: flag for whether to anneal the beta pop parameter

     anneal_beta_compact: flag for whether to anneal the beta compactness parameter

     anneal_beta_segregation: flag for whether to anneal the beta segregation parameter

     anneal_beta_similar: flag for whether to anneal the beta similarity parameter

     adjswap: flag for whether to force algorithm to only do adjacent swaps

     exact_mh: flag for whether to calculate the exact metropolis-hastings w boundary correction

   */

	 // Set nsims if annealing
   int start_anneal;
   int start_cold;
   NumericVector beta_seq(num_annealing_steps);
   if(adapt_beta == "annealing"){
     nsims = num_hot_steps + num_annealing_steps + num_cold_steps;
     start_anneal = num_hot_steps;
     start_cold = num_hot_steps + num_annealing_steps;
     for(int i = 0; i < num_annealing_steps; i++){
       beta_seq(i) = (double)i/num_annealing_steps;
     }
   }

  // Preprocess vector of congressional district assignments
  if(min(cdvec) == 1){
    for(int i = 0; i < cdvec.size(); i++){
      cdvec(i)--;
    }
  }

  // Preprocess vector of original congressional district assignments
  if(min(cdorigvec) == 1){
    for(int i = 0; i < cdorigvec.size(); i++){
      cdorigvec(i)--;
    }
  }

  // Get populations of districts
  NumericVector district_pops = init_pop(popvec, cdvec);

  // Get vector of unique district ids
  NumericVector uniquedists;
  for(int i = 0; i < cdvec.size(); i++){
    if(is_true(any(uniquedists == cdvec(i))) == FALSE){
      uniquedists.push_back(cdvec(i));
    }
  }

  // Get ssd denominator
  double ssd_denom;
  if(beta_compact != 0.0 || anneal_beta_compact == 1){
    ssd_denom = as<double>(calc_betacompact(cdvec,
					    cdvec,
					    popvec,
					    beta_compact,
					    uniquedists,
					    ssdmat,
					    1.0)["compact_new_psi"]);
  }else{
    ssd_denom = 1;
  }

  // Define parity, min and max popoulations
  double parity = sum(popvec) / (max(cdvec) + 1);
  double dist_parity = parity * pct_dist_parity;
  double min_parity = parity - dist_parity;
  double max_parity = parity + dist_parity;

  // Set counter variable
  int k = 0;
  // For storing beta sequence
  int z = 0;
  // For printing progress
  int nsims_10pct = ceil((double)nsims / 10);

  // Store outputted congressional districts
	  if(adapt_beta != "annealing"){
  NumericMatrix cd_store(cdvec.size(), nsims);
  }else{
  NumericMatrix cd_store(cdvec.size(), 1);
	}

  // Store metropolis-hastings decisions for swaps
  NumericVector decision_store(nsims);
  NumericVector mhprob_store(nsims);
  int decision_counter = 0;

  // Store value of psi for all constraints
  NumericVector psipop_store(nsims);
  NumericVector psicompact_store(nsims);
  NumericVector psisegregation_store(nsims);
  NumericVector psisimilar_store(nsims);

  // Store value of p, lambda, weights for all simulations
  NumericVector pparam_store(nsims);

  // Store number of connected components along the boundary, boundary weights
  NumericVector boundarypartitions_store(nsims);
  NumericVector boundaryratio_store(nsims);

  // Store sequence of betas - geyer thompson
  NumericVector betaseq_store(nsims);
  if(anneal_beta_population == 1){
    betaseq_store[z] = beta_population;
  }
  if(anneal_beta_compact == 1){
    betaseq_store[z] = beta_compact;
  }
  if(anneal_beta_segregation == 1){
    betaseq_store[z] = beta_segregation;
  }
  if(anneal_beta_similar == 1){
    betaseq_store[z] = beta_similar;
  }

  // Iterate up z
  z++;

  if(anneal_beta_population == 0 && beta_population != 0.0){
    std::fill(betaseq_store.begin(), betaseq_store.end(), beta_population);
  }
  if(anneal_beta_compact == 0 && beta_compact != 0.0){
    std::fill(betaseq_store.begin(), betaseq_store.end(), beta_compact);
  }
  if(anneal_beta_segregation == 0 && beta_segregation != 0.0){
    std::fill(betaseq_store.begin(), betaseq_store.end(), beta_segregation);
  }
  if(anneal_beta_similar == 0 && beta_similar != 0.0){
    std::fill(betaseq_store.begin(), betaseq_store.end(), beta_similar);
  }

  // Store metropolis-hastings decisions - geyer thompson
  NumericVector decision_betaseq_store(nsims);
  NumericVector mhprob_betaseq_store(nsims);

  // Initialize objects
  List aList_con; NumericVector boundary; List swap_partitions;
  List boundary_partitions; List cutedge_lists; int p; List aList_con_prop;
  NumericVector boundary_prop; List boundary_partitions_prop; int decision;
  List get_constraint; List gt_out; NumericVector cdvec_prop; int i;

	if(adapt_beta == "annealing"){
    Rcout << "---------------------------------" << std::endl;
    Rcout << "-- Simulating at hot temperature." << std::endl;
    Rcout << "---------------------------------" << std::endl;
  }

  // Open the simulations
  while(k < nsims){

    /////////////////////////////////////
    // First: determine boundary cases //
    /////////////////////////////////////
    // Modify aList list for connected components within cd
    aList_con = genAlConn(aList, cdvec);

    // Get vector of boundary units
    boundary = findBoundary(aList, aList_con);

    ///////////////////////////////////////////////////////////////////////////
    // Second: within each congressional district, turn on edges with pr = p //
    ///////////////////////////////////////////////////////////////////////////
    // Continue trying until you get p good swaps
    do{

      // First element is connected adjlist, second element is cut adjlist
      cutedge_lists = cut_edges(aList_con, eprob);

      ////////////////////////////////////////////////////////////////////
      // Third: generate a list of connected components within each cd //
      ///////////////////////////////////////////////////////////////////
      /* List of connected partitions after edgecuts - first element is list of
	 partitions, second element is number of partitions */
      boundary_partitions = bsearch_boundary(cutedge_lists["connectedlist"],
					     boundary);

      ///////////////////////////////////////////////////////////////////////
      // Fourth - select several connected components w/ unif distribution //
      ///////////////////////////////////////////////////////////////////////
      // Draw parameter p (number of swaps for iteration of alg) from pois(lambda)
      p = draw_p(lambda);

      // Loop over p, draw p connected components
      swap_partitions = make_swaps(boundary_partitions["bsearch"],
				   aList,
				   cdvec,
				   cdorigvec,
				   popvec,
				   district_pops,
				   grouppopvec,
				   ssdmat,
				   min_parity,
				   max_parity,
				   p,
				   eprob,
				   beta_population,
				   beta_compact,
				   beta_segregation,
				   beta_similar,
				   ssd_denom);

    }while(as<int>(swap_partitions["goodprop"]) == 0);

    // Get new boundary, then get number of partitions
    if(exact_mh == 1){
      aList_con_prop = genAlConn(aList, as<NumericVector>(swap_partitions["proposed_partition"]));
      boundary_prop = findBoundary(aList, aList_con_prop);
      boundary_partitions_prop = bsearch_boundary(cutedge_lists["connectedlist"],
						  boundary_prop);

      // Correct npartitions to only include boundary partitions that don't break contiguity
      int nvalid_current = count_valid(aList, boundary_partitions["bsearch"], cdvec);
      int nvalid_prop = count_valid(aList, boundary_partitions_prop["bsearch"],
				    swap_partitions["proposed_partition"]);

      // Construct multiple swaps term
      double p_0;
      double F_pi;
      double F_pi_prime;
      if(lambda > 0){
      	p_0 = R::ppois(0, lambda, 1, 0);
      	F_pi = R::ppois(nvalid_current, lambda, 1, 0) - p_0;
      	F_pi_prime = R::ppois(nvalid_prop, lambda, 1, 0) - p_0;
      }else{
      	F_pi = 1.0;
      	F_pi_prime = 1.0;
      }

      // Modify metropolis-hastings ratio
      swap_partitions["mh_prob"] = as<double>(swap_partitions["mh_prob"]) *
	pow((double)nvalid_current / nvalid_prop, (double)p) * (F_pi / F_pi_prime);
      boundaryratio_store(k) = pow((double)nvalid_current / nvalid_prop, (double)p);
    }

    //////////////////////////////////////////
    // Fifth - Accept with some probability //
    //////////////////////////////////////////
    decision = mh_decision(as<double>(swap_partitions["mh_prob"]));

    /////////////////////////////////////////////////////////////
    // Also - for simulated tempering, propose a possible swap //
    /////////////////////////////////////////////////////////////
    if((anneal_beta_population == 1) || (beta_population != 0.0)){ // If constraining, get constraints
      get_constraint = calc_betapop(cdvec,
				    as<NumericVector>(swap_partitions["proposed_partition"]),
				    popvec,
				    beta_population,
				    uniquedists);

      // Store psi value
      if(decision == 1){
	psipop_store[k] = as<double>(get_constraint["pop_new_psi"]);
      }else{
	psipop_store[k] = as<double>(get_constraint["pop_old_psi"]);
      }

      if(anneal_beta_population == 1){ // Tempering step
	if(decision == 1){ // Using value of psi if accepted

	  // Propose swapping beta
	  gt_out = changeBeta(beta_sequence,
			      beta_population,
			      as<double>(get_constraint["pop_new_psi"]),
			      beta_weights,
			      adjswap);

	}else{ // Using value of psi if not accepted

	  // Propose swapping beta
	  gt_out = changeBeta(beta_sequence,
			      beta_population,
			      as<double>(get_constraint["pop_old_psi"]),
			      beta_weights,
			      adjswap);

	}

	// Change beta
	beta_population = as<double>(gt_out["beta"]);

	// Store the output of geyer thompson
	if(k < nsims){
	  betaseq_store[z] = beta_population;
	}
	decision_betaseq_store[k] = as<int>(gt_out["mh_decision"]);
	mhprob_betaseq_store[k] = as<double>(gt_out["mh_prob"]);

      }

    }
    if((anneal_beta_compact == 1) || (beta_compact != 0.0)){ // If constraining, get value of constraints
      get_constraint = calc_betacompact(cdvec,
					as<NumericVector>(swap_partitions["proposed_partition"]),
					popvec,
					beta_compact,
					uniquedists,
					ssdmat,
					ssd_denom);

      // Get psi value
      if(decision == 1){
	psicompact_store[k] = as<double>(get_constraint["compact_new_psi"]);
      }else{
	psicompact_store[k] = as<double>(get_constraint["compact_old_psi"]);
      }

      if(anneal_beta_compact == 1){ // Annealing step
	if(decision == 1){ // Using value of psi if accepted

	  // Propose swapping beta
	  gt_out = changeBeta(beta_sequence,
			      beta_compact,
			      as<double>(get_constraint["compact_new_psi"]),
			      beta_weights,
			      adjswap);

	}else{ // Using value of psi if not accepted

	  // Propose swapping beta
	  gt_out = changeBeta(beta_sequence,
			      beta_population,
			      as<double>(get_constraint["compact_old_psi"]),
			      beta_weights,
			      adjswap);

	}

	// Change beta
	beta_compact = as<double>(gt_out["beta"]);

	// Store the output of geyer thompson
	if(k < nsims){
	  betaseq_store[z] = beta_compact;
	}
	decision_betaseq_store[k] = as<int>(gt_out["mh_decision"]);
	mhprob_betaseq_store[k] = as<double>(gt_out["mh_prob"]);

      }

    }
    if((anneal_beta_segregation == 1) || (beta_segregation != 0.0)){ // If constraining, get value of constraint
      get_constraint = calc_betasegregation(cdvec,
					    as<NumericVector>(swap_partitions["proposed_partition"]),
					    popvec,
					    beta_segregation,
					    uniquedists,
					    grouppopvec);

      // Get psi value
      if(decision == 1){
	psisegregation_store[k] = as<double>(get_constraint["segregation_new_psi"]);
      }else{
	psisegregation_store[k] = as<double>(get_constraint["segregation_old_psi"]);
      }

      if(anneal_beta_segregation == 1){ // Annealing step
	if(decision == 1){

	  // Propose swapping beta
	  gt_out = changeBeta(beta_sequence,
			      beta_segregation,
			      as<double>(get_constraint["segregation_new_psi"]),
			      beta_weights,
			      adjswap);

	}else{ // Use value of psi if not accepted

	  // Propose swapping beta
	  gt_out = changeBeta(beta_sequence,
			      beta_segregation,
			      as<double>(get_constraint["segregation_old_psi"]),
			      beta_weights,
			      adjswap);

	}

	// Change beta
	beta_segregation = as<double>(gt_out["beta"]);

	// Store output of geyer thompson
	if(k < nsims){
	  betaseq_store[k] = beta_segregation;
	}
	decision_betaseq_store[k] = as<int>(gt_out["mh_decision"]);
	mhprob_betaseq_store[k] = as<double>(gt_out["mh_prob"]);

      }

    }
    if((anneal_beta_similar == 1) || (beta_similar != 0.0)){ // If constraining, get value of constraint
      get_constraint = calc_betasimilar(cdvec,
					as<NumericVector>(swap_partitions["proposed_partition"]),
					cdorigvec,
					beta_similar,
					uniquedists);

      // Get psi value
      if(decision == 1){
	psisimilar_store[k] = as<double>(get_constraint["similar_new_psi"]);
      }else{
	psisimilar_store[k] = as<double>(get_constraint["similar_old_psi"]);
      }

      if(anneal_beta_similar == 1){ // Annealing step
	if(decision == 1){

	  // Propose swapping beta
	  gt_out = changeBeta(beta_sequence,
			      beta_similar,
			      as<double>(get_constraint["similar_new_psi"]),
			      beta_weights,
			      adjswap);

	}else{ // Use value of psi if not accepted

	  // Propose swapping beta
	  gt_out = changeBeta(beta_sequence,
			      beta_similar,
			      as<double>(get_constraint["similar_new_psi"]),
			      beta_weights,
			      adjswap);

	}

	// Change beta
	beta_similar = as<double>(gt_out["beta"]);

	// Store output of geyer thompson
	if(k < nsims){
	  betaseq_store[k] = beta_similar;
	}
	decision_betaseq_store[k] = as<int>(gt_out["mh_decision"]);
	mhprob_betaseq_store[k] = as<double>(gt_out["mh_prob"]);

      }

    }
    //////////////////////////////////////
    // Six = clean up and store results //
    //////////////////////////////////////
    cdvec_prop = clone(as<NumericVector>(swap_partitions["proposed_partition"]));
    if(decision == 1){
      // Update cds to proposed cds
      cdvec = clone(as<NumericVector>(swap_partitions["proposed_partition"]));
      // Update district_pops to proposed district pops
      district_pops = clone(as<NumericVector>(swap_partitions["updated_cd_pops"]));
      // Store number of boundary partitions
      boundarypartitions_store[k] = boundary_partitions["npartitions"];
    }else{
      boundarypartitions_store[k] = boundarypartitions_store[k-1];
    }

		// Store previous iteration
    if(adapt_beta != "annealing"){
      for(i = 0; i < cdvec.size(); i++){
					cd_store[k * cdvec.size() + i] = cdvec(i);
      }
    }

    // Store p
    pparam_store[k] = p;

    // Store the decision
    decision_store[k] = decision;
    decision_counter += decision;

    mhprob_store[k] = as<double>(swap_partitions["mh_prob"]);

    // Advance k, z
    k++;
    z++;

    // Print Progress
    if(k % nsims_10pct == 0){
      Rcout << (double)k / nsims_10pct * 10 << " percent done." << std::endl;
      if(adapt_lambda == 1){
	Rcout << "Lambda: " << lambda << std::endl;
      }
      if(adapt_eprob == 1){
	Rcout << "Edgecut Probability: " << eprob << std::endl;
      }
      Rcout << "Metropolis acceptance ratio: "<< (double)decision_counter / (k-1) << std::endl << std::endl;
    }
		if(adapt_beta == "annealing"){
			if(k == start_anneal){
	Rcout << "----------------------------" << std::endl;
	Rcout << "-- Starting annealing stage." << std::endl;
	Rcout << "----------------------------" << std::endl;
			}
			if(k == start_cold){
	Rcout << "----------------------------------" << std::endl;
	Rcout << "-- Simulating at cold temperature." << std::endl;
	Rcout << "----------------------------------" << std::endl;
			}
		}
    // Change eprob, lambda if adaptive
    if(adapt_eprob == 1 || adapt_lambda == 1){
      if(k % 50 == 0){
	if((double)decision_counter / (k-1) > .4){
	  if(adapt_lambda == 1 && lambda < floor((double)aList.size() / 10)){
	    lambda++;
	  }
	  if(adapt_eprob == 1 && eprob < .5){
	    eprob = eprob + .01;
	  }
	}
	if((double)decision_counter / (k-1) < .2){
	  if(adapt_lambda == 1 && lambda > 0){
	    lambda--;
	  }
	  if(adapt_eprob == 1 && eprob > 0){
	    eprob = eprob - .01;
	  }
	}
      }
    }

  }

  // Get distance from parity of each partition
	NumericVector dist_parity_vec;
  NumericVector dist_orig_vec;
  if(adapt_beta != "annealing"){
    dist_parity_vec = distParity(cd_store, popvec);
    dist_orig_vec = diff_origcds(cd_store, cdorigvec);
  }else{
    for(int i = 0; i < cdvec.size(); i++){
      cd_store[i] = cdvec(i);
    }
    dist_parity_vec = distParity(cd_store, popvec);
    dist_orig_vec = diff_origcds(cd_store, popvec);
  }


	// Create list, store output
  List out;
  if(adapt_beta != "annealing"){
    out["partitions"] = cd_store;
    out["distance_parity"] = dist_parity_vec;
    out["distance_original"] = dist_orig_vec;
    out["mhdecisions"] = decision_store;
    out["mhprob"] = mhprob_store;
    out["pparam"] = pparam_store;
    out["beta_sequence"] = betaseq_store;
    out["energy_psi"] = energy_store;
    out["constraint_pop"] = psipop_store;
    out["constraint_compact"] = psicompact_store;
    out["constraint_segregation"] = psisegregation_store;
    out["constraint_similar"] = psisimilar_store;
    out["constraint_countysplit"] = psicountysplit_store;
    out["boundary_partitions"] = boundarypartitions_store;
    out["boundaryratio"] = boundaryratio_store;
    if(adapt_beta == "tempering"){
      out["mhdecisions_beta"] = decision_betaseq_store;
      out["mhprob_beta"] = mhprob_betaseq_store;
    }
    if(adapt_eprob == 1){
      out["final_eprob"] = eprob;
    }
    if(adapt_lambda == 1){
      out["final_lambda"] = lambda;
    }
  }else{
    out["partitions"] = cdvec;
    out["distance_parity"] = dist_parity_vec[0];
    out["distance_original"] = dist_orig_vec[0];
    out["mhdecisions"] = decision_store[k-1];
    out["mhprob"] = mhprob_store[k-1];
    out["pparam"] = pparam_store[k-1];
    out["beta_sequence"] = 1.0;
    out["energy_psi"] = energy_store[k-1];
    out["constraint_pop"] = psipop_store[k-1];
    out["constraint_compact"] = psicompact_store[k-1];
    out["constraint_segregation"] = psisegregation_store[k-1];
    out["constraint_similar"] = psisimilar_store[k-1];
    out["constraint_countysplit"] = psicountysplit_store[k-1];
    out["boundary_partitions"] = boundarypartitions_store[k-1];
    out["boundaryratio"] = boundaryratio_store[k-1];
    if(adapt_eprob == 1){
      out["final_eprob"] = eprob;
    }
    if(adapt_lambda == 1){
      out["final_lambda"] = lambda;
    }
  }

  return out;

}
