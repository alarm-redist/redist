#include "smc_base.h"

// [[Rcpp::export]]
NumericMatrix polsbypopper(IntegerVector from,
                           IntegerVector to,
                           NumericVector area,
                           NumericVector perimeter,
                           IntegerMatrix dm,
                           int nd) {

    NumericMatrix ret(nd, dm.ncol());
    NumericVector dist_area(nd);
    NumericVector dist_peri(nd);
    NumericVector dist_peri2(nd);
    NumericVector zerovec(nd);
    int ne = from.size();
    double pi4 = 4.0*3.14159265;

    for(int c = 0; c < dm.ncol(); c++){
        // update holders:
        dist_area = clone(zerovec);
        dist_peri = clone(zerovec);

        // Get a vector of areas ~ just sum
        for(int r = 0; r < dm.nrow(); r++){
            dist_area(dm(r,c) - 1) += area(r);
        }
        // Get a vector of perims ~ sum by id'ing borders
        for(int e = 0; e < ne; e++){
            if(from(e) == -1){
                dist_peri(dm(to(e) - 1, c) - 1) += perimeter(e);
            } else {
                if(dm(from(e) - 1, c) != dm(to(e) - 1, c)){
                    dist_peri(dm(to(e) - 1, c) - 1) += perimeter(e);
                }
            }
        }

        dist_peri2 = pow(dist_peri, 2.0);
        ret(_,c) = pi4*dist_area/dist_peri2;
    }

    return ret;
}
