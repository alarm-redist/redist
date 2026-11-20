///////////////////////////////////////////////
// Author: Ben Fifield
// Institution: Princeton University
// Date Created: 2014/12/26
// Date Last Modified: 2015/02/26
// Purpose: Contains functions to run calculate beta constraints
///////////////////////////////////////////////

// Header files

#include "redist_types.h"
#include <Rcpp.h>

/* Function to modify adjacency list to reflect adjacency only within
   a particular congressional district */
// [[Rcpp::export]]
Rcpp::List genAlConn(Rcpp::List aList, Rcpp::NumericVector cds) {

    /* Inputs to function:
       aList: adjacency list of geographic units

       cds: vector of congressional district assignments
    */

    // Initialize container list
    Rcpp::List alConnected(cds.size());

    // Initialize
    int i;
    Rcpp::NumericVector avec;
    int cd_i;
    int j;

    // Loop through precincts
    for (i = 0; i < cds.size(); i++) {

        // For precinct i, get adjacent precincts
        avec = aList(i);

        // Get precinct i's congressional district
        cd_i = cds(i);

        // Initialize empty vector
        Rcpp::NumericVector avec_cd;

        // Loop through avec to identify which are in same cd
        for (j = 0; j < avec.size(); j++) {

            // Check if j'th entry in avec is same cd, add to avec_cd if so
            if (cds(avec(j)) == cd_i) {
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
Rcpp::NumericVector findBoundary(Rcpp::List fullList, Rcpp::List conList) {

    /* Inputs to function:
       fullList: Full adjacency list of geographic units

       conList: Adjacency list of geographic units within cong district
    */

    // Initialize container vector of 0's (not boundary) and 1's (boundary)
    Rcpp::NumericVector isBoundary(fullList.size());

    // Initialize inside loop
    Rcpp::NumericVector full;
    Rcpp::NumericVector conn;
    int i;

    // Loop through aList
    for (i = 0; i < fullList.size(); i++) {

        // Get vectors of full and cd-connected components for precinct i
        full = fullList(i);
        conn = conList(i);

        // Compare lengths - if conn < full, then boundary unit
        if (full.size() > conn.size()) {
            isBoundary(i) = 1;
        }
    }

    return isBoundary;
}

