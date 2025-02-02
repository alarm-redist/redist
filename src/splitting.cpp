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
void estimate_cut_k(const SplittingSchedule &splitting_schedule,
    const Graph &g, int &k, int const last_k, 
                      const std::vector<double> &unnormalized_weights, double thresh,
                      double tol, std::vector<std::unique_ptr<Plan>> const &plan_ptrs_vec, 
                      const uvec &counties,
                      Multigraph &cg, const uvec &pop,
                      int const min_region_cut_size, int const max_region_cut_size,
                      bool split_district_only,
                      double const target, int const verbosity) {
    // sample some spanning trees and compute deviances
    int V = g.size();
    int k_max = std::min(10 + (int) (2.0 * V * tol), last_k + 4); // heuristic
    int N_max = plan_ptrs_vec.size();
    int N_adapt = std::min(60 + (int) std::floor(5000.0 / sqrt((double)V)), N_max);

    double lower = target * (1 - tol);
    double upper = target * (1 + tol);

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
    
    for (int i = 0; i < N_max && idx < N_adapt; i++, idx++) {
        if (unnormalized_weights.at(i) == 0) { // skip if not valid
            idx--;
            continue;
        }

        int n_vtx = V;

        // Get the index of the region with the largest dval
        int biggest_region_id; int biggest_dval; int max_valid_dval;

        // if split district only just do remainder 
        if(split_district_only){
            biggest_region_id = plan_ptrs_vec.at(i)->remainder_region;
            biggest_dval = plan_ptrs_vec.at(i)->region_sizes(biggest_region_id);
            max_valid_dval = 1; // max valid dval is 1
        }else{
            biggest_region_id = plan_ptrs_vec.at(i)->region_sizes.head(plan_ptrs_vec.at(i)->num_regions).index_max();
            biggest_dval = plan_ptrs_vec.at(i)->region_sizes(biggest_region_id);
            max_valid_dval = biggest_dval-1;
        }

        int biggest_dval_region_pop = plan_ptrs_vec.at(i)->region_pops.at(biggest_region_id);


        for (int j = 0; j < V; j++) {
            // if not the biggest region mark as ignore
            if (plan_ptrs_vec.at(i)->region_ids(j) != biggest_region_id) {
                ignore[j] = true;
                n_vtx--;
            }
        }

        if (n_vtx > max_V) max_V = n_vtx;

        clear_tree(ust);
        int result = sample_sub_ust(g, ust, V, root, visited, ignore,
                                    pop, lower, upper, counties, cg);
        if (result != 0) {
            idx--;
            continue;
        }

        // reset the cut below pop to zero
        std::fill(cut_below_pop.begin(), cut_below_pop.end(), 0);
        tree_pop(ust, root, pop, cut_below_pop, parents);


        devs.push_back(
        get_ordered_tree_cut_devs(ust, root, cut_below_pop, target, 
                    plan_ptrs_vec.at(i)->region_ids,
                    biggest_region_id, biggest_dval, biggest_dval_region_pop,
                    splitting_schedule.all_regions_min_and_max_possible_cut_sizes[biggest_dval][0],
                    splitting_schedule.all_regions_min_and_max_possible_cut_sizes[biggest_dval][1],
                    splitting_schedule.all_regions_smaller_cut_sizes_to_try[biggest_dval]
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
            int rand_index = r_int(k);
            if(rand_index >= devs.at(i).size()){
                REprintf("For k=%d In est_cut_k rand_index = %d  bigger than devs %d size %d\n",
                k, rand_index, i, (int)  devs.at(i).size());
                throw Rcpp::exception("Erorr in estimate cut k!\n");
            }
            double dev = devs.at(i).at(rand_index);
            if (dev > tol) continue;
            else n_ok++;
            // need min to avoid indexing errors
            for (int j = 0; j < std::min(N_adapt, (int) devs.size()); j++) {
                if(devs.at(j).size() < k) throw Rcpp::exception("Potential k is bigger than region!");
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





