
#include <RcppArmadillo.h>
#include "redist_types.h"
#include "make_swaps_helper.h"
#include "sw_mh_helper.h"
#include "constraint_calc_helper.h"

using namespace Rcpp;

// [[Rcpp::export]]
List rsg(List adj_list,
               NumericVector population,
               int Ndistrict,
               double target_pop,
               double thresh,
               int maxiter
               ) {

	List alConnected;
	int Nprecinct = adj_list.size();
	double maxpop = target_pop * (1+thresh);
	double minpop = target_pop * (1-thresh);
	int i, j, i_dist=0, j_dist=0;
	int p, p_index;
	int iter = 0;
	IntegerVector p_neighbors, d_neighbors;
	IntegerVector idist_pmembers, jdist_pmembers, idist_pneighbors, idist_dneighbors, idist_newmembers;
	NumericVector district_pop(Ndistrict);	//Rcpp Reference says it will always be 0
	NumericVector maxpop_dist(Ndistrict);
	int j_candidates_size, p_neighbors_size;
	IntegerVector idist_members, jdist_members;
	arma::vec i_dist_vec, j_dist_vec, maxpop_dist_vec, j_candidates_vec, p_index_vec;

	// member_dvec tracks which district each precinct belongs to
	// Always of length Nprecinct
	IntegerVector member_dvec(Nprecinct);

	// List of vectors, where each vector is an Integervector of precincts in that district
	// This eventually shrinks to length Ndistrict
	List member_plist(Nprecinct);

	// Step 1: Each precinct represents a single district
	int Ndist_current = Nprecinct;
	for(i=0; i < Nprecinct; i++){
		member_dvec[i] = i;
		IntegerVector onevec(1);
		onevec[0] = i;
		member_plist[i] = onevec;
	}

	while(Ndist_current > Ndistrict){

	  //STEP 2A: Select one of the N (Ndist_current) districts and denote it as district i
	  i_dist_vec = runif(1, 0, 1000000000);
	  i_dist = fmod(i_dist_vec(0), Ndist_current);

	  // Identify neighboring precincts (assume each precinct has at least one neighbor)
	  // Loop over all member precincts to add their adjacency lists
	  // Remove the precincts that are in district(idist) at the end
	  idist_pmembers = IntegerVector(member_plist[i_dist]);
	  idist_pneighbors = IntegerVector(adj_list[idist_pmembers(0)]);
	  if( idist_pmembers.size() > 1 ){
	    for(i=1; i < idist_pmembers.size(); i++){
	      idist_pneighbors = union_( idist_pneighbors, IntegerVector(adj_list[idist_pmembers(i)]) );
	    }
	  } //end if
	  idist_pneighbors = setdiff( idist_pneighbors, idist_pmembers );

	  // Generate list of adjacent districts
	  // Get each neighboring precinct's district membership, and take it's unique vector
	  // Note: union_(x,x) gives the unique elements of x
	  idist_dneighbors = IntegerVector(idist_pneighbors.size());
	  for(i=0; i < idist_dneighbors.size(); i++){
	    idist_dneighbors[i] = member_dvec[ idist_pneighbors[i] ];
	  }
	  idist_dneighbors = unique(idist_dneighbors);

	  //STEP 2B: Select one of district i's bordering districts at random and denote it as district j
	  j_dist_vec = runif(1, 0, 1000000000);
	  j_dist = idist_dneighbors[fmod(j_dist_vec(0), idist_dneighbors.size())];
	  jdist_pmembers = IntegerVector(member_plist[j_dist]);

	  //STEP 2C: Merge district i with district j to form a new district
	  // We merge j into i --- so ensure that i < j and swap if necessary
	  if( i_dist > j_dist){
	    i = i_dist;
	    i_dist = j_dist;
	    j_dist = i;
	  }

	  // Loop over all precincts, update member_dvec 2 ways
	  // Every jdist becomes idist
	  // jdist gets wiped, so decrement every district counter above it
	  for(i=0; i< Nprecinct; i++){
	    if(member_dvec[i] == j_dist) member_dvec[i] = i_dist;
	    if(member_dvec[i] > j_dist) member_dvec[i] = member_dvec[i] - 1;
	  }

	  // Update member_plist
	  idist_newmembers = union_(idist_pmembers, jdist_pmembers) ;
	  member_plist[i_dist] = idist_newmembers;	// merge member lists in idist
	  member_plist.erase(j_dist);				//delete jdist's membership list

	  Ndist_current = Ndist_current - 1;

	}	// end while(Ndist_current > Ndistrict)


	// Calculate district population
	// Only total this once --- we iteratively only add/subtract a single precinct from this after each switch
	for(i=0; i< Nprecinct; i++)	district_pop[ member_dvec[i] ] += population[i];


	iter=0;
	while( (max(district_pop) > maxpop) || (min(district_pop) < minpop) ){

	  // WHILE POPULATION CONSTRAINT NOT MET

	  // STEP 3A: Identify the most populated district and denote it as district i
	  i_dist = which_max(district_pop);

	  // If more than one district has max population, break the tie randomly
	  j=0;
	  for(i=0; i<Ndistrict; i++){
	    if(district_pop[i] == district_pop[i_dist]){
	      maxpop_dist[j] = i;
	      j++;
	    }
	  }
	  maxpop_dist_vec = runif(1, 0, 1000000000);
	  i_dist = maxpop_dist[fmod(maxpop_dist_vec(0), j)];

	  // STEP 3B: Select a random precinct in that district, and denote it as precinct p
	  idist_pmembers = IntegerVector(member_plist[i_dist]);
	  p_index_vec = runif(1, 0, 1000000000);
	  p_index = fmod(p_index_vec(0), idist_pmembers.size());
	  p = idist_pmembers[p_index];


	  // STEP 3C: If precinct p can be reassigned from district i to a new district
	  // without violating contiguity of either this district or district i, reassign it
	  p_neighbors = IntegerVector(adj_list[p]);
	  p_neighbors_size = p_neighbors.size();
	  d_neighbors = IntegerVector(p_neighbors_size);
	  for(i=0; i< p_neighbors_size; i++)	d_neighbors[i] = member_dvec[p_neighbors[i]];


	  // Generate list of adjacent districts
	  IntegerVector j_candidates;
	  j_candidates = unique(d_neighbors);

	  // Decrementing downwards instead to ensure erasing doesn't change full coverage
	  for(i = (j_candidates.size() - 1); i > -1; i--){
		if( j_candidates[i] == i_dist ) j_candidates.erase(i);
	  }
	  j_candidates_size = j_candidates.size();

	  // If there is a valid move, choose a valid j_dist to move to
	  if( j_candidates_size > 0){

	    j_candidates_vec = runif(1, 0, 1000000000);
	    j_dist = j_candidates[fmod(j_candidates_vec(0), j_candidates_size)];

		//Check check that move from i_dist to j_dist is valid
		NumericVector member_new_dvec;
		member_new_dvec = NumericVector(Nprecinct);
		for(i=0; i<Nprecinct; i++) member_new_dvec[i] = member_dvec[i];
	    member_new_dvec[p] = j_dist;
		alConnected = genAlConn(adj_list, member_new_dvec);

		//Make the move if valid
		if(countpartitions(alConnected) == Ndistrict){

			district_pop[i_dist] -= population[p];
	    	district_pop[j_dist] += population[p];

	    	member_dvec[p] = j_dist;

	    	//Assumed here that if idist_members.size==1 (a single precinct district), it will never be the largest district
	    	idist_members = IntegerVector(member_plist[i_dist]);
		    for(i= (idist_members.size() - 1); i > -1 ; i--){
		      if(idist_members[i] == p) idist_members.erase(i);
		    }
		    jdist_members = IntegerVector(member_plist[j_dist]);
		    jdist_members.push_back(p);

		    member_plist[i_dist] = idist_members;
		    member_plist[j_dist] = jdist_members;
			}

	  } // end if(i_valid==1){

	  iter++;
	  if(iter > maxiter){
	    List result;
	    result["plan"] = NumericVector::get_na();
	    result["district_list"] = NumericVector::get_na();
	    result["district_pop"] = NumericVector::get_na();
	    return(result);
	  }

	}  // end while( max(district_pop) > maxpop){

	// Output is:
	// member_dvec - Vector of length #precincts with CD assigments
	// member_plist - List of length #districts with precinct ID's
	// district_pop - vector of length #districts with district populations
	List result;
	result["plan"] = member_dvec;
	result["district_list"] = member_plist;
	result["district_pop"] = district_pop;

	return(result) ;
}


