/////////////////////////////////////
// Author: Ben Fifield
// Institution: Princeton University
// Date Created: 2014/12/17
// Date Modified: 2015/02/26
// Purpose: Supporting functions for swMH() function in redist
/////////////////////////////////////

// Header files
#include "sw_mh_helper.h"

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

  // Initialize
  int i; int pop; arma::uvec cd_i_ind; int j;

  // Loop through cd assignments
  for(i = 0; i < ncds; i++){

    // Initialize population count
    pop = 0;

    // Get indices of cds
    cd_i_ind = find(cds == i);

    // Loop through cd_i_ind, get population values
    for(j = 0; j < cd_i_ind.n_elem; j++){
      pop += popvec(cd_i_ind(j));
    }

    // Put in distpop
    distpop(i) = pop;

  }

  return distpop;

}

// Function to make unidirectional adjacency list bidirectional
List add_ties(List aList){

  // Initialize
  int i; NumericVector list1; int j; NumericVector list2;

  // Loop through vectors in aList
  for(i = 0; i < aList.size(); i++){

    // Get i'th entry in list
    list1 = aList(i);

    // Loop through elements in list1
    for(j = 0; j < list1.size(); j++){

      // Extract adjacency vector for j'th element of i's adjacency list
      list2 = aList(list1(j));

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
  double threshold_prob = 1.0 - eprob;

  // Define lists to store cut-edge and uncut-edge vectors
  List aList_uncut(aList_con.size());
  List aList_cut(aList_con.size());

  // Initialize inside loop
  int i; NumericVector cc_vec_i_all; NumericVector cc_vec_i;
  arma::vec draws;

  // Define list to store output of both lists

  // Loop through elements of aList_con
  for(i = 0; i < aList_con.size(); i++){

    // Extract i'th vector in list
    cc_vec_i_all = aList_con(i);

    // Subset cc_vec_i to elements > i
    cc_vec_i = cc_vec_i_all[cc_vec_i_all > i];

    // For each element in vector, take random draw from [0,1] uniform
    draws = runif(cc_vec_i.size());

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

  // Initialize objects inside loop
  int u; bool in_part; NumericVector adj_u; int i; int v;

  // Begin do{} loop - run until number of elements in boundary_indices is 0
  do{

    // Begin while{} loop - run until q is empty
    while(q.size() > 0){
      Rcpp::checkUserInterrupt();
      // Dequeue first element in queue
      u = q(0);

      // Mark that element in ledger
      mark(u) = u;

      // Check if element is in the partition - add to partition if false
      in_part = is_true(any(partition == u));
      if(in_part == false){
        partition.push_back(u);
      }

      // Get adjacency vector for unit u
      adj_u = aList(u);

      // Loop through elements of adj_u, add to queue and mark if not reached
      if(adj_u.size() > 0){

        // Start loop
        for(i = 0; i < adj_u.size(); i++){

          // Reach element v
          v = adj_u(i);

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
      for(i = boundary_indices.n_elem - 1; i >= 0; i--){
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

  // Get breadth search size
  int bsearch_size = bsearch.size();

  // Get weight_boundary vector
  double weight_boundary = (double)countpartitions(aList) / bsearch_size;

  List out;
  out["bsearch"] = bsearch;
  out["npartitions"] = bsearch_size;
  out["weight_boundary"] = weight_boundary;

  return out;

}

/* Function to count number of valid partitions to swap */
int count_valid(List aList, List boundarypart, NumericVector cdvec){

  int cd_boundary; arma::vec part; int j; int i;
  arma::uvec find_cds; int counter = 0;

  for(i = 0; i < boundarypart.size(); i++){

    // Get the partition
    part = as<arma::vec>(boundarypart(i));

    // Get the congressional district of the boundary
    cd_boundary = cdvec(part(0));

    // Find indices within that congressional district
    find_cds = find(as<arma::vec>(cdvec) == cd_boundary);

    // Remove elements in the partition from that cd
    NumericVector cd_less_boundary;
    for(j = 0; j < find_cds.n_elem; j++){
      if(any(part == find_cds(j)) == false){
        cd_less_boundary.push_back(find_cds(j));
      }
    }

    // If cd_less_boundary empty, then continue
    // Eliminates district so invalid partition
    if(cd_less_boundary.size() == 0){
      continue;
    }

    // Create new adjacency list
    List newadj(cd_less_boundary.size());
    for(j = 0; j < newadj.size(); j++){

      // Extract vector from adjacency list
      NumericVector getadjvec = aList(cd_less_boundary(j));

      // Subset down to elements in cd_less_boundary
      NumericVector getadjvec_sub;
      for(int k = 0; k < getadjvec.size(); k++){
        if(any(as<arma::vec>(cd_less_boundary) == getadjvec(k))){
          getadjvec_sub.push_back(getadjvec(k));
        }
      }

      // Change indices
      NumericVector getadjvec_new;
      for(int k = 0; k < getadjvec_sub.size(); k++){
        arma::uvec ind = find(as<arma::vec>(cd_less_boundary) ==
          getadjvec_sub(k));
        getadjvec_new.push_back(ind(0));
      }

      // Add to newadj
      newadj(j) = getadjvec_new;

    }

    // Calculate number of partitions
    int nparts = countpartitions(newadj);
    if(nparts == 1){
      counter++;
    }

  }

  return counter;

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
                NumericVector pop_vec,
                NumericVector cd_pop_vec,
                List constraints,
                CharacterVector psi_names,
                double minparity,
                double maxparity,
                double parity,
                int p,
                double eprob,
                double beta,
                const Graph &g)
{

  /* Inputs to function:
   boundary_cc: Connected components on district boundaries

   aList: full adjacency list

   cds_old: Current cong district assignments

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

   beta_vra: strength of constraint for segregating subgroup

   beta_partisan: strength of constraint for minimizing partisan bias measure

   beta_similar: strength of constraint for similarity to orig plan

   ssd_denominator: normalizing constant for sum of squared distance psi

   partisan_measure: string, only "efficiency-gap" implemented thus far

   */

  // Initialize objects for swap //
  NumericVector cds_prop = clone(cds_old);
  NumericVector cdspop_prop = clone(cd_pop_vec);
  NumericVector accepted_partitions;
  NumericVector prop_cd_pops;
  NumericVector cds_test;

  // Initialize metropolis-hastings probabilities
  double mh_prob = 1.0;

  NumericVector old_psi(psi_names.size());
  old_psi.names() = psi_names;
  NumericVector new_psi(psi_names.size());
  new_psi.names() = psi_names;

  // Number of unique congressional districts
  int ndists = max(cds_old) + 1;
  int nprec = cds_old.size();

  // Break indicators
  int breakp = 0;
  int numaccept = 0;
  int goodprop = 0;

  IntegerVector curr_cd_swaps(p);
  IntegerVector prop_cd_swaps(p);
  int n_swaps = 0;

  // Begin loop over p
  for(int i = 0; i < p; i++){

    // While there are still possible components available
    int breakwhile = 0;
    int curr_cd;
    int prop_cd;
    NumericVector prop_partitions;
    while(boundary_cc.size() > 0){
      Rcpp::checkUserInterrupt();

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
    curr_cd_swaps(i) = curr_cd;
    prop_cd_swaps(i) = prop_cd;

    // Update metropolis-hastings probabilities
    mh_prob = update_mhprob(prop_partitions,
                            aList,
                            cds_prop,
                            prop_cd,
                            eprob,
                            mh_prob);

    // Update cd assignments and cd populations
    cds_prop = cds_test;
    cdspop_prop = prop_cd_pops;

    // Push back prop_partition to accepted_partitions
    for(int j = 0; j < prop_partitions.size(); j++){
      accepted_partitions.push_back(prop_partitions(j));
    }

    n_swaps++;

    // If we get enough partitions
    if(breakp == 1){
      break;
    }

  } // end loop over p

  // Now calculate target pieces:
  IntegerVector swaps(2 * n_swaps);

  for (int i = 0; i < n_swaps; i++) {
      swaps(i) = curr_cd_swaps(i);
      swaps(i + n_swaps) = prop_cd_swaps(i);
  }

  swaps = sort_unique(swaps);

  // make them 1 idxed
  swaps = swaps + 1;

  std::vector<int> swaps_v = as<std::vector<int>>(swaps);

  mat districts(cds_prop.size(), 2, fill::zeros);
  for(int r = 0; r < nprec; r++){
      districts(r, 0) = cds_prop(r) + 1;
      districts(r, 1) = cds_old(r) + 1;
  }
  arma::umat udistricts = conv_to<umat>::from(districts);
  arma::uvec pops = conv_to<arma::uvec>::from(as<arma::vec>(pop_vec));

  // Multiply mh_prob by constraint values
  double energy_new = calc_gibbs_tgt(udistricts.col(0), ndists, nprec, swaps_v,
                                     new_psi, pops, parity, g, constraints);
  double energy_old = calc_gibbs_tgt(udistricts.col(1), ndists, nprec, swaps_v,
                                     old_psi, pops, parity, g, constraints);

  mh_prob = (double)mh_prob * exp(-1.0 * beta * (energy_new - energy_old));

      // Create returned list
      List out;
      out["proposed_partition"] = cds_prop;
      out["mh_prob"] = mh_prob;
      out["updated_cd_pops"] = cdspop_prop;
      out["goodprop"] = goodprop;
      out["energy_new"] = energy_new;
      out["energy_old"] = energy_old;
      out["new_psi"] = new_psi;
      out["old_psi"] = old_psi;

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
  double mhprobGT = (double)exp(-1 * constraint * (propBeta - beta)) * wj / wi * qji / qij;
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

