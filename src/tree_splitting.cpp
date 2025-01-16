/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2025/1
* Purpose: Functions for finding an edge in a tree to 
remove (splitting the tree)
********************************************************/

#include "tree_splitting.h"

// Given min and max sizes for one of the regions this returns 
// for loop bounds where you only do each pair once
std::pair<int, int> get_potential_region_size_for_loop_bounds(
    const int total_region_size,
    const int min_potential_cut_size, const int max_potential_cut_size
){
    
    // if biggest possible cut size is leq half just return the same bounds
    if(max_potential_cut_size <= total_region_size/2){
        return std::pair<int, int>(min_potential_cut_size, max_potential_cut_size);
    }else if(min_potential_cut_size > total_region_size/2){
        // else if smallest possible size if more than half just subtract 
        // total_region_size and flip 
        return std::pair<int, int>(
            total_region_size - max_potential_cut_size,
            total_region_size - min_potential_cut_size
            );
    }else{ // else must be true that 
    // min_potential_cut_size <= total_region_size/2 < max_potential_cut_size

    return std::pair<int, int>(
            std::min(min_potential_cut_size, total_region_size - max_potential_cut_size),
            total_region_size/2
            );

    }
    
}


/*
 * Calculate the deviations for cutting at every edge in a spanning tree.
//'
//'
//' For each edge it returns the larger of the two deviations associated with 
//' the best region sizes assignment 
 */
std::vector<double> tree_cut_devs(Tree &ust, int root,
                             std::vector<int> const &cut_below_pop, double const target,
                             const arma::subview_col<arma::uword> &region_ids,
                             int const region_id, int const region_size, int const region_pop,
                             int const min_potential_cut_size, int const max_potential_cut_size
                             ) {
    int V = cut_below_pop.size();
    // get for loop bounds that avoid redundant checks 
    std::pair<int, int> loop_bounds = get_potential_region_size_for_loop_bounds(
    region_size,
    min_potential_cut_size, max_potential_cut_size
    );
    // compile a list of candidate edges to cut
    std::vector<double> devs; 
    devs.reserve(V); // reserve V which is overkill but ok
    // REprintf("Only looking for region id %d!\n", region_id);
    // REprintf("Starting at %d and %2.f and ", region_pop, static_cast<double>(region_pop)*2.0);
    for (int i = 0; i < V; i++) {
        // ignore vertices not in the region
        if (i == root || region_ids(i) != region_id) continue;


        // start at total pop since deviance will never be more than total_pop/2
        double smallest_dev = static_cast<double>(region_pop)*2.0;
        
        int below_pop = cut_below_pop.at(i);
        int above_pop = region_pop - below_pop;

        // if one of the populations is zero give large deviance and continue
        if(below_pop == 0 || above_pop == 0){
            devs.push_back(smallest_dev);
            continue;
        } 


        // for each possible d value get the deviation above and below
        for(int cut_region1_size = loop_bounds.first; cut_region1_size <= loop_bounds.second; cut_region1_size++){
            int cut_region2_size = region_size - cut_region1_size;
            double cut_region1_target = target*cut_region1_size; 
            double cut_region2_target = target*cut_region2_size;


            // REprintf("Size 1=%d-Target %.2f, Size 2=%d-Target %.2f", cut_region1_size, cut_region1_target, cut_region2_size, cut_region2_target);
            // find the larger deviation associated with assigning cut_region1_size to cutting below 
            double max_below_dev = std::max(
                std::fabs(below_pop - cut_region1_target) / cut_region1_target,
                std::fabs(above_pop - cut_region2_target) / cut_region2_target
            );
            // find the larger deviation associated with assigning cut_region1_size to cutting above
            double max_above_dev = std::max(
                std::fabs(above_pop - cut_region1_target) / cut_region1_target,
                std::fabs(below_pop - cut_region2_target) / cut_region2_target
            );
            // REprintf("pair dev %.2f, %.2f", max_below_dev, max_above_dev);

            // take this minimum of this d value and all previous ones 
            double min_dev = std::min(
                max_below_dev,
                max_above_dev
            ) ;
            smallest_dev = std::min(smallest_dev, min_dev);
        }
        // REprintf("Max dev is %.3f\n", smallest_dev);
        devs.push_back(smallest_dev);
    }

    std::sort(devs.begin(), devs.end());

    return devs;
}


//' Attempt to pick one of the top k tree edges to split uniformly at random
//'
//' 
//' Attempts to pick one of the top k tree edges to split uniformly at random
//' and if successful returns information on the edge and region sizes 
//' associated with the cut. This function is based on `cut_districts` in `smc.cpp`
//' however the crucial difference is even if a successful cut is found it does not
//' update the plan or the tree.
//'
//'
//' It will only attempt to create regions where the size is between
//' min_potential_d and max_potential_d (inclusive). So the one district
//' split case is `min_potential_d=max_potential_d=1`.
//' 
//' Best edge here is defined as the smallest deviation where the deviation for
//' an edge, size is defined as abs(population - target*size)/target*size. 
//'
//' 
//' @param root The root vertex of the spanning tree
//' @param pop_below The population corresponding to cutting below each vertex. 
//' So `pop_below[v]` is the population associated with the region made by cutting
//' below the vertex `v`
//' @param tree_vertex_parents The parent of each vertex in the tree. A value of -1
//' means the vertex is the root or it is not in the tree.
//' @param k_param The number of best edges we should choose from uniformly
//' at random
//' @param min_potential_cut_size The smallest potential region size to try for a cut. 
//' @param max_potential_cut_size The largest potential region size it will try for a cut. 
//' Setting this to 1 will result in only 1 district splits. 
//' @param region_ids A vector mapping 0 indexed vertices to their region id number
//' @param region_id_to_split The id of the region in the plan object we're attempting to split
//' @param total_region_pop The total population of the region being split 
//' @param total_region_size The size of the region being split 
//' @param lower Acceptable lower bounds on a valid district's population
//' @param upper Acceptable upper bounds on a valid district's population
//' @param target Ideal population of a valid district. This is what deviance is calculated
//' relative to
//'
//' @details No modifications made
//'
//' @return <True, information on the edge cut> if two valid regions were 
//' successfully split, false otherwise
//'
std::pair<bool, EdgeCut> get_naive_top_k_edge(const int root, 
                     const std::vector<int> &pop_below, const std::vector<int> &tree_vertex_parents,
                     const int k_param, 
                     const int min_potential_cut_size, const int max_potential_cut_size,
                     const arma::subview_col<arma::uword> &region_ids,
                     const int region_id_to_split, const int total_region_pop, const int total_region_size,
                     const double lower, const double upper, const double target){
    
    int V = static_cast<int>(region_ids.n_elem);
    // compile a list of: things for each edge in tree
    std::vector<int> candidates; // candidate edges to cut,
    std::vector<double> deviances; // how far from target pop.
    std::vector<bool> is_ok; // whether they meet constraints
    std::vector<int> new_d_val; // the new d_n,k value

    // Reserve at least as much space for top k_param of them
    candidates.reserve(k_param);
    deviances.reserve(k_param);
    is_ok.reserve(k_param);
    new_d_val.reserve(k_param);

    if(region_ids(root) != region_id_to_split){
        REprintf("Root vertex %u is not in region to split %d!\n", region_ids(root), region_id_to_split);
        throw Rcpp::exception("Root vertex is not in region to split!");
    }

    // Now loop over all valid edges to cut
    for (int i = 1; i <= V; i++) { // 1-indexing here
        // Ignore any vertex not in this region or the root vertex as we wont be cutting those
        if (region_ids(i-1) != region_id_to_split || i - 1 == root) continue;

        // Get the population of one side of the partition removing that edge would cause
        int below = pop_below.at(i - 1);
        // Get the population of the other side
        int above = total_region_pop - below;

        // keeps track of the information for the best d_val cut
        double smallest_dev = target * max_potential_cut_size * max_potential_cut_size; // start off with fake maximal value
        bool is_dev2_bigger = false;
        bool cut_is_ok = false;
        int cut_dval;

        // WATCH OUT FOR TYPE

      // if one of the populations is zero immediately reject
        // if(below_pop == 0 || above_pop == 0){
        //     // REprintf("\n");
        //     continue;
        // } 

        // Now try each potential d_nk value from min_potential_d up to max_potential_d
        for(int potential_d = min_potential_cut_size; potential_d <= max_potential_cut_size; potential_d++){

            // dev1 corresponds to region made by cutting below
            double dev1 = std::fabs(below - target * potential_d) / (target * potential_d);
            // dev2 corresponds to region made by cutting above
            double dev2 = std::fabs(above - target * potential_d) / (target * potential_d);

            // if dev1 is smaller of the two and beats the last one then update
            if(dev1 <= dev2 && dev1 < smallest_dev){
                is_dev2_bigger = true;
                smallest_dev = dev1;
                // check if the split induced fits population bounds
                cut_is_ok = lower * potential_d <= below
                            && below <= upper * potential_d
                            && lower * (total_region_size - potential_d) <= above
                            && above <= upper * (total_region_size - potential_d);

                cut_dval = potential_d;
            }else if(dev2 < dev1 && dev2 < smallest_dev){
                is_dev2_bigger = false;
                smallest_dev = dev2;
                // check if the split induced fits population bounds
                cut_is_ok = lower * potential_d <= above
                            && above <= upper * potential_d
                            && lower * (total_region_size - potential_d) <= below
                            && below <= upper * (total_region_size - potential_d);

                cut_dval = potential_d;
            }
        }

        
        if(is_dev2_bigger){
            candidates.push_back(i); // candidate > 0 means we store size of cut below
        } else{
            candidates.push_back(-i); // candidate < 0 means we store size of cut above
        }

        deviances.push_back(smallest_dev);
        is_ok.push_back(cut_is_ok);
        new_d_val.push_back(cut_dval);
    }

    // if less than k_param candidates immediately reject
    if((int) candidates.size() < k_param){
        return std::make_pair(false, EdgeCut());
    }

    
    int idx = r_int(k_param);
    idx = select_k(deviances, idx + 1);
    int cut_vertex = std::fabs(candidates.at(idx)) - 1;


    // reject sample if not ok
    if (!is_ok.at(idx)){
        return std::make_pair(false, EdgeCut());
    }

    // if ok then get the parent of the cut vertex 
    int cut_vertex_parent = tree_vertex_parents.at(cut_vertex);
    int cut_below_pop = pop_below.at(cut_vertex);
    int cut_above_pop = total_region_pop - cut_below_pop;

    int cut_below_region_size; // The size of the region below made by cutting 
    int cut_above_region_size; // The size of the region above made by cutting 

    if (candidates.at(idx) > 0) { // Means we stored size of cut below region
        cut_below_region_size = new_d_val.at(idx);
        cut_above_region_size = total_region_size - cut_below_region_size;
    } else { // Means we stored size of cut above region
        cut_above_region_size = new_d_val.at(idx);
        cut_below_region_size = total_region_size - cut_above_region_size;
    }

    // Now create the EdgeCut object
    EdgeCut selected_edge_cut(root, cut_vertex, cut_vertex_parent,
            cut_below_region_size, cut_below_pop,
            cut_above_region_size, cut_above_pop);

    return std::make_pair(true, selected_edge_cut);

}



//' Pick one of the valid tree edges to split uniformly at random if possible
//'
//' 
//' Returns a valid tree edge to split uniformly at random if at least one
//' valid edge to cut is in the tree. If successful returns information on the 
//' edge and region sizes associated with the cut. 
//'
//' Note even if a successful cut is found it does not
//' update the plan or the tree.
//'
//'
//' It will only attempt to create regions where the size is between
//' min_potential_d and max_potential_d (inclusive). So the one district
//' split case is `min_potential_d=max_potential_d=1`.
//' 
//' Valid edge here is defined as an edge and region sizes such that the 
//' two induced regions both fall within the population bounds.
//'
//' 
//' @param root The root vertex of the spanning tree
//' @param pop_below The population corresponding to cutting below each vertex. 
//' So `pop_below[v]` is the population associated with the region made by cutting
//' below the vertex `v`
//' @param tree_vertex_parents The parent of each vertex in the tree. A value of -1
//' means the vertex is the root or it is not in the tree.
//' @param k_param The number of best edges we should choose from uniformly
//' at random
//' @param min_potential_cut_size The smallest potential region size to try for a cut. 
//' @param max_potential_cut_size The largest potential region size it will try for a cut. 
//' Setting this to 1 will result in only 1 district splits. 
//' @param region_ids A vector mapping 0 indexed vertices to their region id number
//' @param region_id_to_split The id of the region in the plan object we're attempting to split
//' @param total_region_pop The total population of the region being split 
//' @param total_region_size The size of the region being split 
//' @param lower Acceptable lower bounds on a valid district's population
//' @param upper Acceptable upper bounds on a valid district's population
//' @param target Ideal population of a valid district. This is what deviance is calculated
//' relative to
//'
//' @details No modifications made
//'
//' @return <True, information on the edge cut> if two valid regions were 
//' successfully split, false otherwise
//'
std::pair<bool, EdgeCut> get_unif_valid_edge(const int root, 
                     const std::vector<int> &pop_below, const std::vector<int> &tree_vertex_parents,
                     const int min_potential_cut_size, const int max_potential_cut_size,
                     const arma::subview_col<arma::uword> &region_ids,
                     const int region_id_to_split, const int total_region_pop, const int total_region_size,
                     const double lower, const double upper, const double target){
    
    int V = static_cast<int>(region_ids.n_elem);
    // compile a list of: things for each edge in tree
    std::vector<int> candidates; // candidate edges to cut,
    std::vector<int> new_d_val; // the new d_n,k value

    if(region_ids(root) != region_id_to_split){
        REprintf("Root vertex %u is not in region to split %d!\n", region_ids(root), region_id_to_split);
        throw Rcpp::exception("Root vertex is not in region to split!");
    }

    // get for loop bounds that avoid redundant checks 
    std::pair<int, int> loop_bounds = get_potential_region_size_for_loop_bounds(
    total_region_size,
    min_potential_cut_size, max_potential_cut_size
    );

    // Now loop over all valid edges to cut
    for (int i = 1; i <= V; i++) { // 1-indexing here
        // Ignore any vertex not in this region or the root vertex as we wont be cutting those
        if (region_ids(i-1) != region_id_to_split || i - 1 == root) continue;

        // Get the population of one side of the partition removing that edge would cause
        int below = pop_below.at(i - 1);
        // Get the population of the other side
        int above = total_region_pop - below;

        // Now just try every potential below size 
        for(int potential_below_region_size = min_potential_cut_size; 
                potential_below_region_size <= max_potential_cut_size; 
                potential_below_region_size++){
            // check if assigning potential_region_size to cut below leads to valid region
            bool cut_below_ok = lower * potential_below_region_size <= below
                            && below <= upper * potential_below_region_size
                            && lower * (total_region_size - potential_below_region_size) <= above
                            && above <= upper * (total_region_size - potential_below_region_size);

            if(cut_below_ok){
                candidates.push_back(i); // candidate > 0 means we store size of cut below
                new_d_val.push_back(potential_below_region_size);
                }
            }


        // Now try each potential size value from min_potential_d up to max_potential_d
        for(int potential_region_size = loop_bounds.first; 
                potential_region_size <= loop_bounds.second; 
                potential_region_size++){
            
            // check if assigning potential_region_size to cut below leads to valid region
            bool cut_below_ok = lower * potential_region_size <= below
                            && below <= upper * potential_region_size
                            && lower * (total_region_size - potential_region_size) <= above
                            && above <= upper * (total_region_size - potential_region_size);

            if(cut_below_ok){
                candidates.push_back(i); // candidate > 0 means we store size of cut below
                new_d_val.push_back(potential_region_size);
            }

            // check if assigning potential_region_size to cut above leads to valid region
            bool cut_above_ok = lower * potential_region_size <= above
                            && above <= upper * potential_region_size
                            && lower * (total_region_size - potential_region_size) <= below
                            && below <= upper * (total_region_size - potential_region_size);

            // to avoid double counting 
            if(cut_above_ok && (total_region_size - potential_region_size) != potential_region_size){
                candidates.push_back(-i); // candidate < 0 means we store size of cut above
                new_d_val.push_back(potential_region_size);
            }
        }
    }

    // if no valid candidates reject immediately 
    if(candidates.size() == 0){
        return std::make_pair(false, EdgeCut());
    }

    // else select an index uniformly at random
    int idx = r_int(candidates.size());
    int cut_vertex = std::fabs(candidates.at(idx)) - 1;

    // if ok then get the parent of the cut vertex 
    int cut_vertex_parent = tree_vertex_parents.at(cut_vertex);
    int cut_below_pop = pop_below.at(cut_vertex);
    int cut_above_pop = total_region_pop - cut_below_pop;

    int cut_below_region_size; // The size of the region below made by cutting 
    int cut_above_region_size; // The size of the region above made by cutting 

    if (candidates.at(idx) > 0) { // Means we stored size of cut below region
        cut_below_region_size = new_d_val.at(idx);
        cut_above_region_size = total_region_size - cut_below_region_size;
    } else { // Means we stored size of cut above region
        cut_above_region_size = new_d_val.at(idx);
        cut_below_region_size = total_region_size - cut_above_region_size;
    }

    // Now create the EdgeCut object
    EdgeCut selected_edge_cut(root, cut_vertex, cut_vertex_parent,
            cut_below_region_size, cut_below_pop,
            cut_above_region_size, cut_above_pop);

    return std::make_pair(true, selected_edge_cut);

}