#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector contiguity(List adj, IntegerVector group) {
    // unlike in geomander, group takes on a value of 1:n, so we adjust for 0:(n-1)
    IntegerVector grp = group - 1;
    IntegerVector choices = sort_unique(grp);
    IntegerVector conncomp(grp.size());
    IntegerVector group_cc(choices.size());
    IntegerVector currgrp(1);
    IntegerVector temp;
    int cc, s, r;

    for(int i = 0; i < grp.size(); i++){
        if(conncomp(i) == 0){
            group_cc(grp(i)) ++;
            cc = group_cc(grp(i));
            conncomp(i) = cc;

            temp = adj(i);
            std::vector<int> reservoir;
            s = 0;
            for(int j = 0; j < temp.size(); j++){
                if(grp(temp(j)) == grp(i) && conncomp(temp(j)) == 0){
                    reservoir.push_back(temp(j));
                    conncomp(temp(j)) = cc;
                    s++;
                }
            }

            if(s > 0){
                r = 0;
                while(r < s){
                    temp = adj(reservoir[r]);
                    for(int j = 0; j < temp.size(); j++){
                        if(grp(temp(j)) == grp(i) && conncomp(temp(j)) == 0){
                            reservoir.push_back(temp(j));
                            conncomp(temp(j)) = cc;
                            s++;
                        }
                    }
                    r++;
                }
            }
        }
    }

    return conncomp;
}
