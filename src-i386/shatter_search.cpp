#include <Rcpp.h>
using namespace Rcpp;


int shatter_search(List adj_list,
                   int p,
                   int i_dist,
                   IntegerVector member_plist_i,
                   IntegerVector dist_assignment){
  
  // target
  IntegerVector target = member_plist_i;
  IntegerVector adj_p = adj_list[p];
  IntegerVector oneprec = IntegerVector::create(NA_REAL);
  IntegerVector search = rep(oneprec,target.size()); 

  int i = 0, j = p;
  // init the search with the first member of p's neighbors that is in i_dist
  while((j == p) & (i < adj_p.size())){
    if(dist_assignment[adj_p[i]] == i_dist){
      j = adj_p[i];
    }
    i++;
  }

  // If there are no starting points, then the precinct is its own district
  // if so, it should never reach here 
  // (assumption: no two single member districts have largest pop diff)
  if(j == p){
    return(0);
  }
  
  // Now we need to build a new adjacency
  // First we need to keep track of what we've searched
  search[0] = p;
  search[1] = j;
  // Create a progress point to track where we are in the search list
  int progress = 1;
  // similarly track the index for the next open space
  int next_open = 2;
  // We don't need to actually search everything (sometimes we do)
  // so we loop until all links have been found or we are unable to find it all
  // this implementation can be changed, but should be faster for larger cases
  while(is_false(all(in(target, search))) & (progress <  next_open) ){
    adj_p = adj_list[search[progress]];
    for(i = 0; i < adj_p.size(); i++){
      if(dist_assignment[adj_p[i]] == i_dist){
        oneprec = IntegerVector::create(adj_p[i]);
        if(is_false(any(in(oneprec, search)))){
          search[next_open] = adj_p[i];
          next_open++;
        }
      }
    }
    progress++;
  }
  

  // if we searched but didnt find everything
  if(progress >= next_open){
    return 0;
  }
  
  // otherwise we must have left the loop because all entries in search were found in target
  return 1;
  
}


bool can_swap(List adj_list, 
              int p, int j_dist, 
              IntegerVector dist_assignment){
  IntegerVector adj_p = adj_list[p];
  int i;
  for(i = 0; i < adj_p.size(); i++){
    if(dist_assignment[adj_p[i]] == j_dist){
      return TRUE;
    }
  }
  return FALSE;
}
