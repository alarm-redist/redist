/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2025/1
* Purpose: Functions for finding an edge in a tree to 
remove (splitting the tree)
********************************************************/

#include "tree_splitting.h"

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

