//#include <RcppArmadillo.h>
#include <Rcpp.h>

// Originally precinct p is in district i_dist
// Function tests whether all of precinct p's neighbors's (if they are members of district i)
// can still be contiguously connected to district i, if p flips to a different district

int check_contiguity(Rcpp::List adj_list,
					 Rcpp::IntegerVector p_neighbors,
					 int p_neighbors_size,
					 Rcpp::IntegerVector d_neighbors,
					 int i_dist,
					 Rcpp::IntegerVector member_dvec
                     ) {
	
	// pmember_vec tracks which district each precinct belongs to
	//arma::ivec member_dvec(Nprecinct);	

	// List of vectors, where each vector is an Integervector of precincts in that district
	//List member_plist(Ndistrict);	
	
	int i, j;

	// p_neighbors are the precincts adjacent to p
	// other_neighbors are the precincts that are adjacent to the precincts adjacent to p
	Rcpp::IntegerVector other_neighbors;
	
//	for(i=0; i< p_neighbors.size(); i++){
	for(i=0; i< p_neighbors_size; i++){

		// Consider only adjacent precincts that are also in district i
//		if( member_dvec[p_neighbors[i]] == i_dist){
		if( d_neighbors[i] == i_dist){

			other_neighbors = adj_list[p_neighbors[i]];
			
			// Looping over neighbors of p_neighbors
			for(j=0; j < other_neighbors.size(); j++){

				if(member_dvec[other_neighbors[j]] == i_dist){
					other_neighbors[j] = 1;
				}
				else{
					other_neighbors[j] = 0;
				}

				// If there is no adjacent precinct other than p that is in district i,
				// contiguity is broken, so reject move
				if( sum(other_neighbors) == 1) return(0);
			
			}  // end for(j=0; j < other_neighbors.size(); j++)

		} // end if( member_dvec[p_neighbors[i]] == i_dist)

	} //end for(i=0; i< p_neighbors_size; i++)

	// No adjacent precinct has it's contiguity broken
	return(1);
}
