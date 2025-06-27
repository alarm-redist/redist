/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2024/11
* Purpose: Functions for Splitting Trees and Plans
********************************************************/

#include "splitting.h"



/*
 * Choose k and multiplier for efficient, accurate sampling
 */
void estimate_cut_k(
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    RNGState &rng_state,
    int &k, int const last_k, 
    const arma::vec &unnormalized_weights, double thresh,
    double tol, std::vector<std::unique_ptr<Plan>> const &plan_ptrs_vec, 
    bool split_district_only,
    int const verbosity) {
    // sample some spanning trees and compute deviances
    int V = map_params.V;
    int k_max = std::min(10 + (int) (2.0 * V * tol), last_k + 4); // heuristic
    int N_max = plan_ptrs_vec.size();
    int N_adapt = std::min(60 + (int) std::floor(5000.0 / sqrt((double)V)), N_max);

    double target = map_params.target;
    double lower = target * (1 - tol);
    double upper = target * (1 + tol);
    int num_regions = plan_ptrs_vec[0]->num_regions;

    std::vector<std::vector<double>> devs;
    devs.reserve(N_adapt);
    vec distr_ok(k_max+1, fill::zeros);
    int root;
    int max_ok = 0;
    std::vector<bool> ignore(V);
    std::vector<bool> visited(V);
    std::vector<int> cut_below_pop(V,0);
    std::vector<int> parents(V);

    int idx = 0;
    int max_V = 0;
    Tree ust = init_tree(V);

    bool any_size_split = splitting_schedule.schedule_type == SplittingSizeScheduleType::AnyValidSizeSMD;
    
    for (int i = 0; i < N_max && idx < N_adapt; i++, idx++) {
        if (unnormalized_weights.at(i) == 0) { // skip if not valid
            idx--;
            continue;
        }

        int n_vtx = V;

        // Get the index of the region with the largest dval
        int biggest_region_id; int biggest_region_size; 

        // if split district only just do remainder 
        if(split_district_only){
            biggest_region_id = plan_ptrs_vec.at(i)->num_regions-1;
            biggest_region_size = plan_ptrs_vec.at(i)->region_sizes[biggest_region_id];
        }else if(any_size_split){
            // else get 
            auto max_it = std::max_element(
                plan_ptrs_vec.at(i)->region_sizes.begin(),
                plan_ptrs_vec.at(i)->region_sizes.begin() + num_regions
            );
            biggest_region_id = std::distance(
                plan_ptrs_vec.at(i)->region_sizes.begin(), 
                max_it
            );
            biggest_region_size = plan_ptrs_vec.at(i)->region_sizes[biggest_region_id];
        }else{// custom size 
            biggest_region_size = -1; biggest_region_id = num_regions -1;
            for (int j = 0; j < num_regions; j++)
            {
                int region_size = plan_ptrs_vec[i]->region_sizes[j];
                // check if valid and bigger
                if(splitting_schedule.valid_region_sizes_to_split[region_size] && 
                   region_size > biggest_region_size){
                    biggest_region_id = j;
                    biggest_region_size = region_size;
                }
            }
        }

        

        int biggest_size_region_pop = plan_ptrs_vec[i]->region_pops[biggest_region_id];
        std::pair<int, int> min_and_max_possible_cut_sizes = splitting_schedule.all_regions_min_and_max_possible_cut_sizes[biggest_region_size];
        int min_possible_cut_size = min_and_max_possible_cut_sizes.first;
        int max_possible_cut_size = min_and_max_possible_cut_sizes.second;
        
        for (int j = 0; j < V; j++) {
            // if not the biggest region mark as ignore
            if (plan_ptrs_vec.at(i)->region_ids[j] != biggest_region_id) {
                ignore[j] = true;
                n_vtx--;
            }
        }
        

        // Rprintf("Tree on region %d of size %d has %d vertices and pop %d!\n",
        // biggest_region_id,
        //     biggest_region_size, n_vtx, plan_ptrs_vec.at(i)->region_pops.at(biggest_region_id));

        // plan_ptrs_vec.at(i)->Rprint();
        // Rprintf("\n\n");

        if (n_vtx > max_V) max_V = n_vtx;

        clear_tree(ust);
        int result = sample_sub_ust(map_params.g, ust, V, root, visited, ignore,
                                    map_params.pop, 
                                    min_possible_cut_size*lower, max_possible_cut_size*upper, 
                                    map_params.counties, map_params.cg, rng_state);
        if (result != 0) {
            idx--;
            continue;
        }

        // reset the cut below pop to zero
        std::fill(cut_below_pop.begin(), cut_below_pop.end(), 0);
        tree_pop(ust, root, map_params.pop, cut_below_pop, parents);

        
        devs.push_back(
        get_ordered_tree_cut_devs(ust, root, cut_below_pop, target, 
                    plan_ptrs_vec.at(i)->region_ids,
                    biggest_region_id, biggest_region_size, 
                    biggest_size_region_pop,
                    min_possible_cut_size, max_possible_cut_size,
                    splitting_schedule.all_regions_smaller_cut_sizes_to_try[biggest_region_size]
                    )
                      );

        int n_ok = 0;
        for(const auto &a_dev : devs.at(idx)){
            if (a_dev <= tol) { // sorted
                n_ok++;
            } else {
                break;
            }
        }

        if (n_ok <= k_max)
            distr_ok(n_ok) += 1.0 / N_adapt;
        if (n_ok > max_ok && n_ok < k_max){
            max_ok = n_ok;
        }
            
        Rcpp::checkUserInterrupt();
    }

    if (idx < N_adapt) N_adapt = idx; // if rejected too many in last step
    // For each k, compute pr(selected edge within top k),
    // among maps where valid edge was selected
    for (k = 1; k <= k_max; k++) {
        // if k > max_v
        double sum_within = 0;
        int n_ok = 0;
        for (int i = 0; i < N_adapt; i++) {
            int rand_index = rng_state.r_int(k);
            if(rand_index >= devs.at(i).size()){
                REprintf("For k=%d In est_cut_k rand_index = %d  bigger than devs %d size %d\n",
                k, rand_index, i, (int)  devs.at(i).size());
                continue;
                // throw Rcpp::exception("Erorr in estimate cut k!\n");
            }
            double dev = devs.at(i).at(rand_index);
            if (dev > tol) continue;
            else n_ok++;
            // need min to avoid indexing errors
            for (int j = 0; j < std::min(N_adapt, (int) devs.size()); j++) {
                if(devs.at(j).size() < k){
                    REprintf("For k=%d its  bigger than the number of devs %d for %d\n",
                k, (int) devs.at(j).size(), i);
                    continue;
                    //throw Rcpp::exception("Potential k is bigger than region!");
                } 
                sum_within += ((double) (dev <= devs.at(j).at(k-1))) / N_adapt;
            }
        }
        if (sum_within / n_ok >= thresh) break;
    }

    if (k >= k_max) {
        if (verbosity >= 3) {
            Rcout << " [maximum hit; falling back to naive k estimator]";
        }
        k = max_ok;
    }

    if (last_k < k_max && k < last_k * 0.6) k = (int) (0.5*k + 0.5*last_k);

    k = std::min(std::max(max_ok + 1, k) + 1 - (distr_ok(k) > 0.99) + (thresh == 1),
                 max_V - 1);
}


