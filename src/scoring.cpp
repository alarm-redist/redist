/********************************************************
 * Author: Philip O'Sullivan
 * Institution: Harvard University
 * Date Created: 2025/02
 * Purpose: Implement constraint and scoring function
 ********************************************************/

 #include "scoring.h"

constexpr bool DEBUG_SCORING_VERBOSE = false;


// helpers


std::vector<std::vector<int>> build_admin_vertex_lists(
    const Graph &g, const arma::uvec &admin_units
){
    // we assume admin units is 1 indexed and if `k` units then values are in `1:k`
    int const num_counties = arma::max(admin_units);
    std::vector<std::vector<int>> admin_vertex_lists(num_counties);
    // nothing if only 1 county 
    if(num_counties == 1) return admin_vertex_lists;
    // else walk through the graph and add each vertex to list for each unit
    int const V = g.size();

    for (int v = 0; v < V; v++)
    {
        int v_admin_unit = admin_units[v]-1;
        admin_vertex_lists[v_admin_unit].push_back(v);
    }
    return admin_vertex_lists;
}

// counts how many administrative unit are split by a plan
// We say a unit is split if there is more than one region inside it
// The region_reindex_vec lets us map different region ids to the same 
// index for the purpose of this function
int count_admin_splits(
    std::vector<std::vector<int>> const &admin_vertex_lists, 
    PlanVector const &region_ids, std::vector<int> const &region_reindex_vec
){
    int split_units = 0;
    int const num_admin_units = admin_vertex_lists.size();

    for (size_t admin_unit_i = 0; admin_unit_i < num_admin_units; admin_unit_i++)
    {
        // Find what region the first vertex in the unit is 
        auto const first_unit_vertex = admin_vertex_lists[admin_unit_i][0];
        auto const first_unit_region = region_reindex_vec[region_ids[first_unit_vertex]];
        // We see if all the other vertices are in the same region 
        for (auto const v: admin_vertex_lists[admin_unit_i]){
            auto const v_unit_region = region_reindex_vec[region_ids[v]];
            if(first_unit_region != v_unit_region){
                ++split_units;
                break;
            }
        }
        
    }

    return split_units;
}



int count_total_admin_splits(
    std::vector<std::vector<int>> const &admin_vertex_lists, 
    std::vector<std::set<int>> &admin_unit_regions,
    PlanVector const &region_ids, std::vector<int> const &region_reindex_vec
){
    int total_splits = 0;
    int const num_admin_units = admin_vertex_lists.size();

    // clear each units sets 
    for (size_t i = 0; i < num_admin_units; i++)
    {
        admin_unit_regions[i].clear();
    }
    

    for (size_t admin_unit_i = 0; admin_unit_i < num_admin_units; admin_unit_i++)
    {
        // First count how many different regions are in that admin unit
        for (auto const v: admin_vertex_lists[admin_unit_i]){
            auto const v_unit_region = region_reindex_vec[region_ids[v]];
            admin_unit_regions[admin_unit_i].insert(v_unit_region);
        }
        // The number of splits is the number of regions in the admin unit minus 1
        total_splits += admin_unit_regions[admin_unit_i].size() - 1;
    }

    return total_splits;
}

/*
 * Given a contiguous administrative units this builds a forest on the 
 * admin units where each tree is a spanning tree on an admin unit along
 with a vector of roots. Lets you traverse all counties in O(V) time and
 space
 */
// Spanning forest on the counties, ie each tree is a tree on a specific county
// roots of each county tree, so [i] is root of tree on county[i+1]
std::pair<Tree,std::vector<int>> build_admin_forest(
    const Graph &g, const arma::uvec &admin_units
){
    // we assume admin units is 1 indexed and if `k` units then values are in `1:k`
    int const num_counties = arma::max(admin_units);
    // nothing if only 1 county 
    if(num_counties == 1) return make_pair(Tree(0), std::vector<int>());
    int const V = g.size();
    // make vector tracking which vertices we've visitied
    std::vector<bool> visited(V, false);
    std::vector<bool> admin_unit_visited(num_counties, false);

    // make a vector tracking the roots of each county tree
    std::vector<int> admin_forest_roots(num_counties);
    // init the forest
    Tree admin_forest(V);

    // Go through graph and build tree on each admin unit
    for (int v = 0; v < V; v++)
    {
        // COUNTIES ARE 1 INDEXED!!
        int v_admin_unit = admin_units[v]-1; 
        // skip if we've visitied this county before 
        if(admin_unit_visited[v_admin_unit]){
            // sanity check can delete later
            if(!visited[v]) throw Rcpp::exception("Should have visitied this vertex already!\n");
            continue;
        } 
        

        // else we can start to build a tree on it 
        // mark this as the root
        admin_forest_roots[v_admin_unit] = v;
        visited[v] = true;
        // queue with vertex 
        std::queue<int> vertex_queue;
        // add the root 
        vertex_queue.push(v);

        // keep going through the children until queue not empty
        while(!vertex_queue.empty()){
            // get from queue
            int u = vertex_queue.front();
            vertex_queue.pop();
            // mark as visited since it has to share this county
            visited[u] = true;
            int u_admin_unit = admin_units[u]-1;
            // sanity check delete later
            if(u_admin_unit != v_admin_unit){
                REprintf("v county %d, u county %d", admin_units(u)-1, v_admin_unit);
                throw Rcpp::exception("County forest went wrong!!\n");
            } 

            // see if any children in the same county
            for (int const child_vertex : g[u]){
                // add to queue if same county and not visitied yet
                if(admin_units[child_vertex]-1 == u_admin_unit && !visited[child_vertex]){
                    // add in tree
                    admin_forest[u].push_back(child_vertex);
                    // mark as visited to avoid being added later  
                    visited[child_vertex] = true;
                    vertex_queue.push(child_vertex);
                }
            }
        }
        // mark this county as visited 
        admin_unit_visited[v_admin_unit] = true;
    }

    return std::make_pair(admin_forest, admin_forest_roots);
}


int count_min_threshold_regions(
    int const num_populations,
    std::vector<arma::vec> const &group_pop, 
    std::vector<arma::vec> const &total_pop, 
    std::vector<double> const &min_fracs,
    std::vector<bool> const &region_ids_to_count,
    PlanVector const &region_ids, std::vector<int> const &region_reindex_vec,
    std::vector<std::vector<double>> &plan_group_pops, 
    std::vector<std::vector<double>> &plan_total_pops
){
    // clear the group populations 
    for (size_t i = 0; i < num_populations; i++)
    {
        std::fill(plan_group_pops[i].begin(), plan_group_pops[i].end(), 0.0);
        std::fill(plan_total_pops[i].begin(), plan_total_pops[i].end(), 0.0);
    }
    


    auto const V = region_ids.size();
    for (size_t v = 0; v < V; v++)
    {
        // get the reinexed region id
        auto const v_region = region_reindex_vec[region_ids[v]];
        // REprintf("Count=%s| v=%u region %u reindex %u\n", (region_ids_to_count[v_region] ? "TRUE" : "FALSE"),
        //     v, region_ids[v], region_reindex_vec[region_ids[v]]);
        // check if this region is one we care about 
        if(!region_ids_to_count[v_region]) continue;

        for (size_t i = 0; i < num_populations; i++)
        {
            // add to the counts for these two regions 
            plan_group_pops[i][v_region] += group_pop[i][v];
            plan_total_pops[i][v_region] += total_pop[i][v];
        }
        
    }

    // Now count how many clear the threshold 
    int regions_above_threshold = 0;

    for (size_t region_id = 0; region_id < region_ids_to_count.size(); region_id++)
    {
        // skip if we don't count this region 
        if(!region_ids_to_count[region_id]) continue; 

        bool region_ok = true;
        
        for (size_t i = 0; i < num_populations; i++)
        {
            // break if population is zero or ratio over threshold
            if(plan_total_pops[i][region_id] == 0.0 ||
                plan_group_pops[i][region_id] / plan_total_pops[i][region_id] < min_fracs[i]){
                region_ok = false;
                break;
            }
        }

        // skip if all regions are ok
        if(region_ok){
            regions_above_threshold++;
        }

        // REprintf("Region %u - Ratio %f\n", region_id, pop_ratio);
    }
    
    return regions_above_threshold;

}


//

std::pair<bool, double> RegionConstraint::compute_region_score(const Plan &plan, int region_id) const{
    double region_score = strength * compute_raw_region_constraint_score(
        plan.num_regions, 
        plan.region_ids, plan.region_sizes, plan.region_pops, 
        region_id);
    
    if(hard_constraint){
        if(region_score >= hard_threshold){
            return std::make_pair(false, region_score);
        }else{
            return std::make_pair(true, 0);
        }
    }else{
        return std::make_pair(true, region_score);
    }
}

double RegionConstraint::compute_soft_region_score(
    int const num_regions, 
    PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
    int const region_id
) const{
    if(hard_constraint){
        return 0.0;
    }else{
        return strength * compute_raw_region_constraint_score(
        num_regions, 
        region_ids, region_sizes, region_pops, 
        region_id);
    }
}


std::pair<bool, double> RegionConstraint::compute_merged_region_score(
    const Plan &plan, int region1_id, int region2_id) const{
    double region_score = strength * compute_raw_merged_region_constraint_score(
        plan.num_regions, 
        plan.region_ids, plan.region_sizes, plan.region_pops, 
        region1_id, region2_id);

    if(hard_constraint){
        if(region_score >= hard_threshold){
            return std::make_pair(false, region_score);
        }else{
            return std::make_pair(true, 0);
        }
    }else{
        return std::make_pair(true, region_score);
    }
}



bool RegionConstraint::region_constraint_ok(const Plan &plan, int region_id) const{
    if (!hard_constraint){
        return true;
    }else{
        double region_score = strength * compute_raw_region_constraint_score(
            plan.num_regions, 
            plan.region_ids, plan.region_sizes, plan.region_pops, 
            region_id
        );
        if(DEBUG_SCORING_VERBOSE){
            REprintf("Score %f, thresh %f so %d\n", region_score, hard_threshold, region_score < hard_threshold);
        }
        return region_score < hard_threshold;
    }
};


 /********************************************************
 * Constraints Supported
 *   - Population tempering (PopTemperConstraint)
 *      This one enoucrages smaller deviations 
 *   - Group Hinge (GroupHingeConstraint)
 *   - Incumbency (IncumbentConstraint)
 *  
 ********************************************************/



double PopTemperConstraint::compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const{
    int const region_pop = region_pops[region_id];
    int const region_size = region_sizes[region_id];
    return compute_log_pop_temper(target, pop_temper, ndists, region_pop, region_size);
}


double PopTemperConstraint::compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const{
    int const merged_region_pop = region_pops[region1_id] + region_pops[region2_id];
    int const merged_region_size = region_sizes[region1_id] + region_sizes[region2_id];
    return compute_log_pop_temper(target, pop_temper, ndists, merged_region_pop, merged_region_size);
}



double PopDevConstraint::compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const{
    double raw_score = eval_pop_dev(
        region_ids, 
        region_id, region_id,
        total_pop, parity
    );

    return raw_score;
}
// log constraint for region made by merging region 1 and 2
double PopDevConstraint::compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const{
    double raw_score = eval_pop_dev(
        region_ids, 
        region1_id, region2_id,
        total_pop, parity
    );
    return raw_score;
}


double StatusQuoConstraint::compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const{
    double raw_score = eval_sq_entropy(
        region_ids, current,
        region_id, region_id,
        pop, 
        ndists, n_current, V
    );

    return raw_score;
}
// log constraint for region made by merging region 1 and 2
double StatusQuoConstraint::compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const{
    double raw_score = eval_sq_entropy(
        region_ids, current,
        region1_id, region2_id,
        pop, 
        ndists, n_current, V
    );
    return raw_score;
}

double SegregationConstraint::compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const{
    double raw_score = eval_segregation(
        region_ids, region_id, region_id,
        V, grp_pop, total_pop
    );

    return raw_score;
}
// log constraint for region made by merging region 1 and 2
double SegregationConstraint::compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const{
    double raw_score = eval_segregation(
        region_ids, region1_id, region2_id,
        V, grp_pop, total_pop
    );
    return raw_score;
}


double GroupPowerConstraint::compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const{

    double raw_score = eval_grp_pow(
        region_ids, V, region_id, region_id,
        grp_pop, total_pop,
        tgt_grp, tgt_other, pow
    );

    return raw_score;
}
// log constraint for region made by merging region 1 and 2
double GroupPowerConstraint::compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const{
    double raw_score = eval_grp_pow(
        region_ids, V, region1_id, region2_id,
        grp_pop, total_pop,
        tgt_grp, tgt_other, pow
    );
    return raw_score;
}

double GroupHingeConstraint::compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const{
    double raw_score = eval_grp_hinge(
        region_ids, V, region_id, region_id,
        tgts_group, group_pop, total_pop
    );

    return raw_score;
}
// log constraint for region made by merging region 1 and 2
double GroupHingeConstraint::compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const{
    double raw_score = eval_grp_hinge(
        region_ids, V, region1_id, region2_id,
        tgts_group, group_pop, total_pop
    );
    return raw_score;
}



double IncumbentConstraint::compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const{
    double raw_score = eval_inc(
        region_ids, 
        region_id, region_id,
        incumbents);

    return raw_score;
}
// log constraint for region made by merging region 1 and 2
double IncumbentConstraint::compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const{
    double raw_score = eval_inc(
        region_ids, 
        region1_id, region2_id,
        incumbents);
    return raw_score;
}





double SplitsConstraint::compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const{
    double raw_score = eval_splits(
        region_ids, region_id,
        admin_units, n_admin_units, smc
    );

    
    return raw_score;
}
// log constraint for region made by merging region 1 and 2
double SplitsConstraint::compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const{
    // make a copy of the plans 
    std::vector<RegionID> dummy_merged_vec(
        region_ids.begin(),
        region_ids.end()
    );
    // now make all instance of region2 region1
    for (size_t i = 0; i < dummy_merged_vec.size(); i++)
    {
        if(dummy_merged_vec[i] == region2_id){
            dummy_merged_vec[i] = region1_id;
        }
    }
    // very crude, if a bottleneck can probably just make eval splits support
    // implicit merging 
    double raw_score = eval_splits(
        dummy_merged_vec, region1_id,
        admin_units, n_admin_units, smc
    );
    return raw_score;
}

double MultisplitsConstraint::compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const{
    double raw_score = eval_multisplits(
        region_ids, region_id,
        admin_units, n_admin_units, smc
    );

    
    return raw_score;
}
// log constraint for region made by merging region 1 and 2
double MultisplitsConstraint::compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const{
    // make a copy of the plans 
    std::vector<RegionID> dummy_merged_vec(
        region_ids.begin(),
        region_ids.end()
    );
    // now make all instance of region2 region1
    for (size_t i = 0; i < dummy_merged_vec.size(); i++)
    {
        if(dummy_merged_vec[i] == region2_id){
            dummy_merged_vec[i] = region1_id;
        }
    }
    // very crude, if a bottleneck can probably just make eval splits support
    // implicit merging 
    double raw_score = eval_multisplits(
        dummy_merged_vec, region1_id,
        admin_units, n_admin_units, smc
    );
    return raw_score;
}


double TotalSplitsConstraint::compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const{
    double raw_score = eval_total_splits(
        region_ids, region_id,
        admin_units, n_admin_units, smc
    );

    return raw_score;
}
// log constraint for region made by merging region 1 and 2
double TotalSplitsConstraint::compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const{
    // make a copy of the plans 
    std::vector<RegionID> dummy_merged_vec(
        region_ids.begin(),
        region_ids.end()
    );
    // now make all instance of region2 region1
    for (size_t i = 0; i < dummy_merged_vec.size(); i++)
    {
        if(dummy_merged_vec[i] == region2_id){
            dummy_merged_vec[i] = region1_id;
        }
    }
    // very crude, if a bottleneck can probably just make eval splits support
    // implicit merging 
    double raw_score = eval_total_splits(
        dummy_merged_vec, region1_id,
        admin_units, n_admin_units, smc
    );
    return raw_score;
}


double PolsbyConstraint::compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const{

    double raw_score = eval_polsby(
        region_ids, region_id, region_id, V,
        from, to, area, perimeter
    );

    return raw_score;
}
// log constraint for region made by merging region 1 and 2
double PolsbyConstraint::compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const{
    double raw_score = eval_polsby(
        region_ids, region1_id, region2_id, V,
        from, to, area, perimeter
    );

    return raw_score;
}


double CustomRegionConstraint::compute_raw_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int region_id) const{
    // Need to copy into Rcpp vector since no SEXP for current region ids
     std::copy(
        region_ids.begin(),
        region_ids.end(),
        rcpp_plan_wrap.begin()
    );

    
    double raw_score = static_cast<double>(
        as<NumericVector>(fn(rcpp_plan_wrap, region_id))[0]
    );

    return raw_score;
}
// log constraint for region made by merging region 1 and 2
double CustomRegionConstraint::compute_raw_merged_region_constraint_score(
            int const num_regions, 
            PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
            int const region1_id, int const region2_id) const{
    // Need to copy into Rcpp vector since no SEXP for current region ids
    // now make all instance of region2 region1
    for (size_t i = 0; i < region_ids.size(); i++)
    {
        if(region_ids[i] == region2_id){
            rcpp_plan_wrap[i] = region1_id;
        }else{
            rcpp_plan_wrap[i] = region_ids[i];
        }
    }
    
    double raw_score = static_cast<double>(
        as<NumericVector>(fn(rcpp_plan_wrap, region1_id))[0]
    );

    return raw_score;
}


void PlanConstraint::print() const{
    REprintf("Scoring Plans with Following Number of Regions:");
    bool first_element = true;
    // just skip the first element bc size is ndists+1
    for (int i = 1; i < num_regions_to_score.size(); i++)
    {
        if(num_regions_to_score[i]) REprintf("%d,", i);
    }

    REprintf("\n");
}


bool PlanConstraint::plan_constraint_ok(const Plan &plan) const{
    if (!hard_constraint){
        return true;
    }else{
        double region_score = strength * compute_raw_plan_constraint_score(
            plan.num_regions, plan.region_ids, plan.region_sizes, plan.region_pops
        );
        if(DEBUG_SCORING_VERBOSE){
        REprintf("Score %f, thresh %f\n", region_score, hard_threshold);
        }
        return region_score < hard_threshold;
    }
}


std::pair<bool, double> PlanConstraint::compute_plan_score(
    int const num_regions, 
    PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops
) const{
    double plan_score = strength * compute_raw_plan_constraint_score(num_regions, region_ids, region_sizes, region_pops);

    if(hard_constraint){
        if(plan_score >= hard_threshold){
            return std::make_pair(false, plan_score);
        }else{
            return std::make_pair(true, 0);
        }
    }else{
        return std::make_pair(true, plan_score);
    }
}


std::pair<bool, double> PlanConstraint::compute_merged_plan_score(
    const Plan &plan, 
    int const region1_id, int const region2_id) 
const{
    double plan_score = strength * compute_raw_merged_plan_constraint_score(plan, region1_id, region2_id);
    if(hard_constraint){
        if(plan_score >= hard_threshold){
            return std::make_pair(false, plan_score);
        }else{
            return std::make_pair(true, 0);
        }
    }else{
        return std::make_pair(true, plan_score);
    }
}

double CustomPlanConstraint::compute_raw_plan_constraint_score(
    int const num_regions, 
    PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops
) const{
    // Need to copy into Rcpp vector since no SEXP for current region ids
    Rcpp::IntegerVector rcpp_plan_wrap(
        region_ids.begin(),
        region_ids.end()
    );

    Rcpp::IntegerVector rcpp_sizes_wrap(
        region_sizes.begin(),
        region_sizes.end()
    );
    
    double raw_score = static_cast<double>(
        as<NumericVector>(fn(rcpp_plan_wrap, rcpp_sizes_wrap, num_regions))[0]
    );

    return raw_score;
}
// log constraint for region made by merging region 1 and 2
double CustomPlanConstraint::compute_raw_merged_plan_constraint_score(const Plan &plan, int const region1_id, int const region2_id) const{
    // Need to copy into Rcpp vector since no SEXP for current region ids
    Rcpp::IntegerVector rcpp_plan_wrap(
        plan.region_ids.begin(),
        plan.region_ids.end()
    );

    Rcpp::IntegerVector rcpp_sizes_wrap(
        plan.region_sizes.begin(),
        plan.region_sizes.end()
    );

    // now make all instance of region2 region1
    for (size_t i = 0; i < plan.region_ids.size(); i++)
    {
        if(rcpp_plan_wrap[i] == region2_id){
            rcpp_plan_wrap[i] = region1_id;
        }
    }

    // adjust the sizes as well so region 1 has all the sizes
    rcpp_sizes_wrap[region1_id] += rcpp_sizes_wrap[region2_id];
    rcpp_sizes_wrap[region2_id] = 0;
    
    double raw_score = static_cast<double>(
        as<NumericVector>(fn(rcpp_plan_wrap, rcpp_sizes_wrap, plan.num_regions-1))[0]
    );

    return raw_score;
}


double PlanSplitsConstraint::compute_raw_plan_constraint_score(
    int const num_regions, 
    PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops
) const{
    // no splits if blank map or only 1 admin unit 
    if(num_regions == 1 || num_admin_units == 1) return 0;

    // set the reindex for each region to be itself
    std::iota(region_reindex_vec.begin(), region_reindex_vec.end(), 0);

    auto splits = count_admin_splits(
        admin_vertex_lists, 
        region_ids, region_reindex_vec
    );

    return splits;

}

double PlanSplitsConstraint::compute_raw_merged_plan_constraint_score(
    const Plan &plan, int const region1_id, int const region2_id) const{
    // no splits if blank map or only 1 admin unit 
    if(plan.num_regions == 2 || num_admin_units == 1) return 0;

    for (int region_id = 0; region_id < plan.num_regions; region_id++)
    {
        if(region_id == region2_id){
            region_reindex_vec[region2_id] = region1_id;
        }else{
            region_reindex_vec[region_id] = region_id;
        }
    }

    auto splits = count_admin_splits(
        admin_vertex_lists, 
        plan.region_ids, region_reindex_vec
    );

    return splits;
}


double TotalPlanSplitsConstraint::compute_raw_plan_constraint_score(
    int const num_regions, 
    PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops
) const{
    // no splits if blank map or only 1 admin unit 
    if(num_regions == 1 || num_admin_units == 1) return 0;

    // set the reindex for each region to be itself
    std::iota(region_reindex_vec.begin(), region_reindex_vec.end(), 0);

    auto splits = count_total_admin_splits(
        admin_vertex_lists, admin_unit_regions,
        region_ids, region_reindex_vec
    );

    return splits;

}

double TotalPlanSplitsConstraint::compute_raw_merged_plan_constraint_score(
    const Plan &plan, int const region1_id, int const region2_id) const{
    // no splits if blank map or only 1 admin unit 
    if(plan.num_regions == 2 || num_admin_units == 1) return 0;

    for (int region_id = 0; region_id < plan.num_regions; region_id++)
    {
        if(region_id == region2_id){
            region_reindex_vec[region2_id] = region1_id;
        }else{
            region_reindex_vec[region_id] = region_id;
        }
    }

    auto splits = count_total_admin_splits(
        admin_vertex_lists, admin_unit_regions,
        plan.region_ids, region_reindex_vec
    );

    return splits;
}


double MinGroupFracConstraint::compute_raw_plan_constraint_score(
    int const num_regions, 
    PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops
) const{
    // if just 1 region we just do the sums
    if(num_regions == 1){
        double pops_above = 0.0;
        for (size_t i = 0; i < num_populations; i++)
        {
            if(arma::sum(group_pops[i])/arma::sum(total_pops[i]) >= min_fracs[i]){
                pops_above++;
            }
        }
        
        return pops_above;
    }

    // set the reindex for each region to be itself
    std::iota(region_reindex_vec.begin(), region_reindex_vec.end(), 0);
    std::fill(region_ids_to_count.begin(), region_ids_to_count.end(), false);

    // set the reindex for each region to be itself and make all of those ok
    for (size_t i = 0; i < num_regions; i++)
    {
        region_reindex_vec[i] = i;
        region_ids_to_count[i] = true;
    }

    
    
    auto num_regions_above_threshold = count_min_threshold_regions(
        num_populations, group_pops, total_pops, 
        min_fracs,
        region_ids_to_count,
        region_ids, region_reindex_vec,
        plan_group_pops, plan_total_pops
    );

    // if(num_regions_above_threshold >= 1){
    //     REprintf("Min Frac %f | %d regions above!\n", min_frac, num_regions_above_threshold);
    // }

    return static_cast<double>(num_regions_above_threshold);

}

double MinGroupFracConstraint::compute_raw_merged_plan_constraint_score(
    const Plan &plan, int const region1_id, int const region2_id) const{

    // if just 1 region we just sum the two 
    if(plan.num_regions == 1){
        throw Rcpp::exception("Calling MinGroupFracConstraint Merge on a 1 region plan!\n");
        return 0.0;
    }else if(plan.num_regions == 2){
        double pops_above = 0.0;
        for (size_t i = 0; i < num_populations; i++)
        {
            if(arma::sum(group_pops[i])/arma::sum(total_pops[i]) >= min_fracs[i]){
                pops_above++;
            }
        }
        
        return pops_above;
    }

    // set the reindex for each region to be itself
    std::iota(region_reindex_vec.begin(), region_reindex_vec.end(), 0);
    std::fill(region_ids_to_count.begin(), region_ids_to_count.end(), false);

    // set the reindex for each region to be itself and make all of those ok
    for (size_t i = 0; i < plan.num_regions; i++)
    {
        region_reindex_vec[i] = i;
        region_ids_to_count[i] = true;
    }

    // merge region2 into region 1 
    region_reindex_vec[region2_id] = region1_id;
    region_ids_to_count[region2_id] = false;
    
    auto num_regions_above_threshold = count_min_threshold_regions(
        num_populations, group_pops, total_pops, 
        min_fracs,
        region_ids_to_count,
        plan.region_ids, region_reindex_vec,
        plan_group_pops, plan_total_pops
    );

    return static_cast<double>(num_regions_above_threshold);
}


double ValidDistrictsConstraint::compute_raw_plan_constraint_score(
    int const num_regions, 
    PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops
) const{
    // REprintf("Valid sizes are:");
    // for (int i = 1; i <= map_params.total_seats; i++)
    // {
    //     if(map_params.is_district[i]) REprintf("%d, ", i);
    // }
    
    // threshold is always .5 so returning 1 means reject
    // first check no region is smaller than the smallest district size 
    for (size_t i = 0; i < num_regions; i++)
    {
        auto region_size = region_sizes[i];
        if(region_size < map_params.smallest_district_size){
            return 1.0;
        }
        // if the plan is all districts, ie ndists = num regions then check every region is a district
        if(num_regions == map_params.ndists && !map_params.is_district[region_size]){
            // if not a district return false
            return 1.0;
        }
    }
    // Now return ok 
    return 0.0;
}

// Scoring function
ScoringFunction::ScoringFunction(
    MapParams const &map_params,
    Rcpp::List const &constraints, double const pop_temper, bool const smc,
    int const thread_id
):
map_params(map_params), num_non_final_soft_region_constraints(0), num_final_soft_region_constraints(0), all_rounds_soft_region_constraints(0), 
total_soft_plan_constraints(0),
num_hard_plan_constraints(0), 
total_soft_constraints(0), num_hard_region_constraints(0),
any_soft_custom_constraints(false), any_hard_custom_constraints(false){
    // First add region constraints 
    // add pop temper if doing that 
    if(pop_temper != 0){
        non_final_region_constraint_ptrs.emplace_back(
            std::make_unique<PopTemperConstraint>(
                map_params.target, map_params.ndists, 
                map_params, pop_temper
            ));
    }

    // Add region constraints 
    // add constraints if in the list 
    if (constraints.containsElementNamed("pop_dev")) {
        Rcpp::List constr = constraints["pop_dev"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            region_constraint_ptrs.emplace_back(
                std::make_unique<PopDevConstraint>(
                    constr_inst, 
                    map_params.target, map_params.pop
                ));
        }
    }
    if (constraints.containsElementNamed("status_quo")) {
        Rcpp::List constr = constraints["status_quo"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            region_constraint_ptrs.emplace_back(
                std::make_unique<StatusQuoConstraint>(
                    constr_inst, 
                    as<arma::uvec>(constr_inst["current"]), map_params.pop,
                    map_params.ndists, as<int>(constr_inst["n_current"]), map_params.V
                ));
        }
    }
    if (constraints.containsElementNamed("segregation")) {
        Rcpp::List constr = constraints["segregation"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            region_constraint_ptrs.emplace_back(
                std::make_unique<SegregationConstraint>(
                    constr_inst, 
                    as<arma::uvec>(constr_inst["group_pop"]), 
                    as<arma::uvec>(constr_inst["total_pop"]), 
                    map_params.V
                ));
        }
    }
    if (constraints.containsElementNamed("grp_pow")) {
        Rcpp::List constr = constraints["grp_pow"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            region_constraint_ptrs.emplace_back(
                std::make_unique<GroupPowerConstraint>(
                    constr_inst, map_params.V,
                    as<arma::uvec>(constr_inst["group_pop"]), 
                    as<arma::uvec>(constr_inst["total_pop"]),
                    as<double>(constr_inst["tgt_group"]),
                    as<double>(constr_inst["tgt_other"]),
                    as<double>(constr_inst["pow"])
                ));
        }
    }
    if (constraints.containsElementNamed("compet")) {
        Rcpp::List constr = constraints["compet"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            // Competition is just group power with group target and other target .5
            arma::uvec dvote = constr_inst["dvote"];
            arma::uvec total = dvote + as<arma::uvec>(constr_inst["rvote"]);
            region_constraint_ptrs.emplace_back(
                std::make_unique<GroupPowerConstraint>(
                    constr_inst, map_params.V,
                    dvote, total,
                    .5, .5, 
                    as<double>(constr_inst["pow"])
            ));
        }
    }
    if (constraints.containsElementNamed("grp_hinge")) {
        // create constraint objects
        Rcpp::List constr = constraints["grp_hinge"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            region_constraint_ptrs.emplace_back(
                std::make_unique<GroupHingeConstraint>(
                    constr_inst, map_params.V, as<arma::vec>(constr_inst["tgts_group"]),
                    as<arma::uvec>(constr_inst["group_pop"]), as<arma::uvec>(constr_inst["total_pop"])
            ));
        }
    }
    if (constraints.containsElementNamed("grp_inv_hinge")) {
        // NOTE: As far as I can tell grp_inv_hinge is same code as grp_inv_hinge
        // create constraint objects
        Rcpp::List constr = constraints["grp_inv_hinge"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            region_constraint_ptrs.emplace_back(
                std::make_unique<GroupHingeConstraint>(
                    constr_inst, map_params.V, as<arma::vec>(constr_inst["tgts_group"]),
                    as<arma::uvec>(constr_inst["group_pop"]), as<arma::uvec>(constr_inst["total_pop"])
                ));
        }
    }
    if (constraints.containsElementNamed("incumbency")) {
        Rcpp::List constr = constraints["incumbency"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            region_constraint_ptrs.emplace_back(
                std::make_unique<IncumbentConstraint>(
                    constr_inst, as<arma::uvec>(constr_inst["incumbents"])
                ));
        }
    }
    if (constraints.containsElementNamed("splits")) {
        Rcpp::List constr = constraints["splits"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            region_constraint_ptrs.emplace_back(
                std::make_unique<SplitsConstraint>(
                    constr_inst, 
                    as<arma::uvec>(constr_inst["admin"]), as<int>(constr_inst["n"]),
                    smc
                ));
        }
    }
    if (constraints.containsElementNamed("multisplits")) {
        Rcpp::List constr = constraints["multisplits"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            region_constraint_ptrs.emplace_back(
                std::make_unique<MultisplitsConstraint>(
                    constr_inst, 
                    as<arma::uvec>(constr_inst["admin"]), as<int>(constr_inst["n"]),
                    smc
                ));
        }
    }
    if (constraints.containsElementNamed("total_splits")) {
        Rcpp::List constr = constraints["total_splits"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            region_constraint_ptrs.emplace_back(
                std::make_unique<TotalSplitsConstraint>(
                    constr_inst, 
                    as<arma::uvec>(constr_inst["admin"]), as<int>(constr_inst["n"]),
                    smc
            ));
        }
    }
    if (constraints.containsElementNamed("polsby")) {
        Rcpp::List constr = constraints["polsby"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            region_constraint_ptrs.emplace_back(
                std::make_unique<PolsbyConstraint>(
                    constr_inst, map_params.V,
                    as<arma::ivec>(constr_inst["from"]), 
                    as<arma::ivec>(constr_inst["to"]),
                    as<arma::vec>(constr_inst["area"]),
                    as<arma::vec>(constr_inst["perimeter"])
                ));
        }
    }

// lp[i] += add_constraint("fry_hold", constraints,
//     [&] (List l) -> double {
//         return eval_fry_hold(districts.col(i), j,
//                              as<uvec>(l["total_pop"]),
//                              as<mat>(l["ssdmat"]),
//                              as<double>(l["denominator"]));
//     });

// lp[i] += add_constraint("qps", constraints,
//     [&] (List l) -> double {
//         return eval_qps(districts.col(i), j,
//                         as<uvec>(l["total_pop"]),
//                         as<uvec>(l["cities"]), as<int>(l["n_city"]),
//                         n_distr);
//     });


    if (constraints.containsElementNamed("custom")) {
        Rcpp::List constr = constraints["custom"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            bool hard_constraint = false;
            if (constr_inst.containsElementNamed("hard_constraint")){
                hard_constraint = as<bool>(constr_inst["hard_constraint"]);
            }
            region_constraint_ptrs.emplace_back(
                std::make_unique<CustomRegionConstraint>(
                    constr_inst, map_params.V,
                    as<Rcpp::Function>(constr_inst["fn"])
                )
            );
            // mark custom R constraints as true 
            any_soft_custom_constraints = true; 
            // if hard custom constraint note that
            if(hard_constraint){
                any_hard_custom_constraints = true;
            }
        }
    }

    // Now add plan constraints 
    // counts splits in the entire plan
    if (constraints.containsElementNamed("plan_splits")) {
        Rcpp::List constr = constraints["plan_splits"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            double strength = constr_inst["strength"];
            std::vector<bool> num_regions_to_score(map_params.ndists+1, true);
            if (constr_inst.containsElementNamed("nregions_to_score")){
                // The vector in R is one indexed but c++ is 0 indexed so need to pad an 
                // extra element
                num_regions_to_score = Rcpp::as<std::vector<bool>>(constr_inst["nregions_to_score"]);
                num_regions_to_score.insert(num_regions_to_score.begin(), false);
            }
            bool hard_constraint = false;
            if (constr_inst.containsElementNamed("hard_constraint")){
                hard_constraint = as<bool>(constr_inst["hard_constraint"]);
            }
            double hard_threshold = 0.0;
            if (constr_inst.containsElementNamed("hard_threshold")){
                hard_threshold = as<double>(constr_inst["hard_threshold"]);
            }
            if (strength != 0) {
                // build the forest and get the roots 
                arma::uvec admin_units = as<arma::uvec>(constr_inst["admin"]);
                auto admin_vertex_lists = build_admin_vertex_lists(map_params.g, admin_units);
                plan_constraint_ptrs.emplace_back(
                    std::make_unique<PlanSplitsConstraint>(
                        strength, map_params.ndists,
                        admin_units, admin_vertex_lists,
                        num_regions_to_score,
                        hard_constraint, hard_threshold
                    )
                );
            }
        }
    }
    // total splits for whole plan
    if (constraints.containsElementNamed("total_plan_splits")) {
        Rcpp::List constr = constraints["total_plan_splits"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            double strength = constr_inst["strength"];
            std::vector<bool> num_regions_to_score(map_params.ndists+1, true);
            if (constr_inst.containsElementNamed("nregions_to_score")){
                // The vector in R is one indexed but c++ is 0 indexed so need to pad an 
                // extra element
                num_regions_to_score = Rcpp::as<std::vector<bool>>(constr_inst["nregions_to_score"]);
                num_regions_to_score.insert(num_regions_to_score.begin(), false);
            }
            bool hard_constraint = false;
            if (constr_inst.containsElementNamed("hard_constraint")){
                hard_constraint = as<bool>(constr_inst["hard_constraint"]);
            }
            double hard_threshold = 0.0;
            if (constr_inst.containsElementNamed("hard_threshold")){
                hard_threshold = as<double>(constr_inst["hard_threshold"]);
            }
            if (strength != 0) {
                // build the forest and get the roots 
                arma::uvec admin_units = as<arma::uvec>(constr_inst["admin"]);
                auto admin_vertex_lists = build_admin_vertex_lists(map_params.g, admin_units);
                plan_constraint_ptrs.emplace_back(
                    std::make_unique<TotalPlanSplitsConstraint>(
                        strength, map_params.ndists,
                        admin_units, admin_vertex_lists,
                        num_regions_to_score,
                        hard_constraint, hard_threshold
                    )
                );
            }
        }
    }
    // counts the number of regions greater than or equal to a certain fraction
    if (constraints.containsElementNamed("min_group_frac")) {
        Rcpp::List constr = constraints["min_group_frac"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            double strength = constr_inst["strength"];
            std::vector<bool> num_regions_to_score(map_params.ndists+1, true);
            if (constr_inst.containsElementNamed("nregions_to_score")){
                // The vector in R is one indexed but c++ is 0 indexed so need to pad an 
                // extra element
                num_regions_to_score = Rcpp::as<std::vector<bool>>(constr_inst["nregions_to_score"]);
                num_regions_to_score.insert(num_regions_to_score.begin(), false);
            }
            bool hard_constraint = false;
            if (constr_inst.containsElementNamed("hard_constraint")){
                hard_constraint = as<bool>(constr_inst["hard_constraint"]);
            }
            double hard_threshold = 0.0;
            if (constr_inst.containsElementNamed("hard_threshold")){
                hard_threshold = as<double>(constr_inst["hard_threshold"]);
            }
            if (strength != 0) {
                plan_constraint_ptrs.emplace_back(
                    std::make_unique<MinGroupFracConstraint>(
                        strength, map_params.ndists, map_params.is_district,
                        as<std::vector<arma::vec>>(constr_inst["group_pops"]),
                        as<std::vector<arma::vec>>(constr_inst["total_pops"]), 
                        as<std::vector<double>>(constr_inst["min_fracs"]),
                        as<double>(constr_inst["num_populations"]),
                        num_regions_to_score,
                        hard_constraint, hard_threshold
                    )
                );
            }

        }
    }

    // Add custom plan constraints 
    if (constraints.containsElementNamed("custom_plan")) {
        Rcpp::List constr = constraints["custom_plan"];
        for (int i = 0; i < constr.size(); i++) {
            List constr_inst = constr[i];
            double strength = constr_inst["strength"];
            std::vector<bool> num_regions_to_score(map_params.ndists+1, true);
            if (constr_inst.containsElementNamed("nregions_to_score")){
                // The vector in R is one indexed but c++ is 0 indexed so need to pad an 
                // extra element
                num_regions_to_score = Rcpp::as<std::vector<bool>>(constr_inst["nregions_to_score"]);
                num_regions_to_score.insert(num_regions_to_score.begin(), false);
            }
            bool hard_constraint = false;
            if (constr_inst.containsElementNamed("hard_constraint")){
                hard_constraint = as<bool>(constr_inst["hard_constraint"]);
            }
            double hard_threshold = 0.0;
            if (constr_inst.containsElementNamed("hard_threshold")){
                hard_threshold = as<double>(constr_inst["hard_threshold"]);
            }
            if (strength != 0) {
                plan_constraint_ptrs.emplace_back(
                    std::make_unique<CustomPlanConstraint>(
                        strength, 
                        as<Rcpp::Function>(constr_inst["fn"]),
                        num_regions_to_score,
                        hard_constraint, hard_threshold
                    )
                );
                any_soft_custom_constraints = true; 
                // if hard custom constraint note that
                if(hard_constraint){
                    any_hard_custom_constraints = true;
                }
            }

        }
    }

    // Now add hard constraints 
    if (constraints.containsElementNamed("plan_valid_district_sizes")) {
        // Add this check 
        plan_constraint_ptrs.emplace_back(
            std::make_unique<ValidDistrictsConstraint>(map_params)
        );
    }


    for(auto const &constraint_ptr: region_constraint_ptrs){
        if(constraint_ptr->hard_constraint){
            ++num_hard_region_constraints;
        }else{
            ++all_rounds_soft_region_constraints;
        }
    }
    for(auto const &constraint_ptr: non_final_region_constraint_ptrs){
        if(constraint_ptr->hard_constraint){
            ++num_hard_region_constraints;
        }else{
            ++num_non_final_soft_region_constraints;
        }
    }

    for(auto const &constraint_ptr: plan_constraint_ptrs){
        if(constraint_ptr->hard_constraint){
            ++num_hard_plan_constraints;
        }else{
            ++total_soft_plan_constraints;
        }
    }

    total_soft_region_constraints = num_non_final_soft_region_constraints+num_final_soft_region_constraints+all_rounds_soft_region_constraints;

    total_soft_constraints = total_soft_region_constraints + total_soft_plan_constraints;

    any_soft_region_constraints = total_soft_region_constraints > 0;
    any_soft_plan_constraints = total_soft_plan_constraints > 0;

    any_hard_plan_constraints = num_hard_plan_constraints > 0;
    any_hard_region_constraints = num_hard_region_constraints > 0;
    any_hard_constraints = static_cast<bool>(any_hard_plan_constraints + any_hard_region_constraints);

    if(DEBUG_SCORING_VERBOSE){
        REprintf(
            "Constraint has %d soft region constraints, %d soft plan ones\n",
            total_soft_region_constraints, total_soft_plan_constraints
        );
    }

}



std::pair<bool, double> ScoringFunction::compute_region_full_score(
    const Plan &plan, 
    int const region_id, bool const is_final
) const{
    // start out with score of zero
    double region_score = 0.0;

    // check if its a multidistrict 
    bool const is_multidistrict = !map_params.is_district[plan.region_sizes[region_id]];

    if(DEBUG_SCORING_VERBOSE){
        REprintf("Scoring region %d, size %u! %s DISTRICT, %s DISTRICT \n", 
            region_id, plan.region_sizes[region_id], 
            (is_multidistrict ? "NOT A" : "IS A"),
            (map_params.is_district[plan.region_sizes[region_id]] ? "IS A" : "NOT A" ));
    }

    // add the score from each constraint 
    for(auto const &constraint_ptr: region_constraint_ptrs){
        // skip if multidistrict and we only score districts
        if(constraint_ptr->score_districts_only && is_multidistrict) continue;

        auto score_result = constraint_ptr->compute_region_score(plan, region_id);
        // immediately return if threshold triggered
        if(!score_result.first){
            return std::make_pair(false, 0.0);
        }else{
            if(DEBUG_SCORING_VERBOSE){
            REprintf("Adding %f!\n", score_result.second);
            }
            region_score += score_result.second;
        }        
    }
    // if final then return 
    if(is_final) return std::make_pair(true, region_score);

    // if not final then add those constraints (for now just pop_temper)
    for(auto const &constraint_ptr: non_final_region_constraint_ptrs){
        // skip if multidistrict and we only score districts
        if(constraint_ptr->score_districts_only && is_multidistrict) continue;

        auto score_result = constraint_ptr->compute_region_score(plan, region_id);
        // immediately return if threshold triggered
        if(!score_result.first){
            return std::make_pair(false, 0.0);
        }else{
            region_score += score_result.second;
        }        
    }

    return std::make_pair(true, region_score);
}


double ScoringFunction::compute_region_soft_score(const Plan &plan, int const region_id, bool const is_final) 
    const{    
    
    return compute_region_full_score(plan, region_id, is_final).second;

}


std::pair<bool, double> ScoringFunction::compute_merged_region_full_score(const Plan &plan, 
    int const region1_id, int const region2_id, bool const is_final) 
    const{
    // start out with log score of zero
    double region_score = 0.0;

    // check if its a multidistrict 
    bool const is_multidistrict = !map_params.is_district[
        plan.region_sizes[region1_id] + plan.region_sizes[region2_id]
    ];

    // add the score from each constraint 
    for(auto const &constraint_ptr: region_constraint_ptrs){
        // skip if multidistrict and we only score districts
        if(constraint_ptr->score_districts_only && is_multidistrict) continue;

        auto score_result = constraint_ptr->compute_merged_region_score(plan, region1_id, region2_id);
        // immediately return if threshold triggered
        if(!score_result.first){
            return std::make_pair(false, 0.0);
        }else{
            region_score += score_result.second;
        }        
    }
    // if final then return 
    if(is_final) return std::make_pair(true, region_score);

    // if not final then add those constraints (for now just pop_temper)
    for(auto const &constraint_ptr: non_final_region_constraint_ptrs){
        // skip if multidistrict and we only score districts
        if(constraint_ptr->score_districts_only && is_multidistrict) continue;

        auto score_result = constraint_ptr->compute_merged_region_score(plan, region1_id, region2_id);
        // immediately return if threshold triggered
        if(!score_result.first){
            return std::make_pair(false, 0.0);
        }else{
            region_score += score_result.second;
        }        
    }

    return std::make_pair(true, region_score);
}


std::pair<bool, double> ScoringFunction::compute_plan_score(const Plan &plan) const{
    double plan_score = 0.0;

    auto const num_regions = plan.num_regions;

    // add the score from each constraint 
    for(auto const &constraint_ptr: plan_constraint_ptrs){
        if(!constraint_ptr->num_regions_to_score[num_regions]) continue;

        auto score_result = constraint_ptr->compute_plan_score(
            plan.num_regions, plan.region_ids, plan.region_sizes, plan.region_pops
        );
        // immediately return if threshold triggered
        if(!score_result.first){
            return std::make_pair(false, 0.0);
        }else{
            plan_score += score_result.second;
        }

    }

    return std::make_pair(true, plan_score);
}


std::pair<bool, double> ScoringFunction::compute_merged_plan_score(
    const Plan &plan, 
    int const region1_id, int const region2_id, bool const is_final) 
const{
    double plan_score = 0.0;
    auto const merged_num_regions = plan.num_regions - 1;

    // add the score from each constraint 
    for(auto const &constraint_ptr: plan_constraint_ptrs){
        if(!constraint_ptr->num_regions_to_score[merged_num_regions]) continue; 

        auto score_result = constraint_ptr->compute_merged_plan_score(plan, region1_id, region2_id);
        // immediately return if threshold triggered
        if(!score_result.first){
            return std::make_pair(false, 0.0);
        }else{
            plan_score += score_result.second;
        }
    }

    return std::make_pair(true, plan_score);
}

bool ScoringFunction::merged_region_ok(
    Plan const &plan, 
    int const region1_id, int const region2_id, 
    bool const is_final_split
) const{
    if(!any_hard_region_constraints) return true;

    // check if its a multidistrict 
    bool const is_multidistrict = !map_params.is_district[
        plan.region_sizes[region1_id] + plan.region_sizes[region2_id]
    ];

    // add the score from each constraint 
    for(auto const &constraint_ptr: region_constraint_ptrs){
        // skip if multidistrict and we only score districts
        if(constraint_ptr->score_districts_only && is_multidistrict) continue;

        auto score_result = constraint_ptr->compute_merged_region_score(plan, region1_id, region2_id);
        // immediately return if threshold triggered
        if(!score_result.first){
            return false;
        }      
    }
    // if final then return 
    if(is_final_split) return true;

    // if not final then add those constraints (for now just pop_temper)
    for(auto const &constraint_ptr: non_final_region_constraint_ptrs){
        // skip if multidistrict and we only score districts
        if(constraint_ptr->score_districts_only && is_multidistrict) continue;

        auto score_result = constraint_ptr->compute_merged_region_score(plan, region1_id, region2_id);
        // immediately return if threshold triggered
        if(!score_result.first){
            return false;
        }            
    }

    return true;

}


bool ScoringFunction::entire_merged_plan_constraint_only_ok(
    Plan const &plan, 
    int const region1_id, int const region2_id, 
    bool const is_final_split
) const{
    if(!any_hard_plan_constraints) return true;

    // this will never be true for a merged plan
    auto const merged_num_regions = plan.num_regions - 1;

    // check each hard constraint 
    for(auto const &constraint_ptr: plan_constraint_ptrs){
        if(!constraint_ptr->hard_constraint || !constraint_ptr->num_regions_to_score[merged_num_regions]) continue; 

        auto score_result = constraint_ptr->compute_merged_plan_score(plan, region1_id, region2_id);
        // immediately return if threshold triggered
        if(!score_result.first){
            return false;
        }
    }

    return true;

}

bool ScoringFunction::merged_plan_ok(Plan const &plan, 
    int const region1_id, int const region2_id, 
    bool const is_final_split) const{

    if(!any_hard_constraints) return true;
    
    if(!merged_region_ok(plan, region1_id, region2_id, is_final_split)){
        return false;
    }

    if(!entire_merged_plan_constraint_only_ok(plan, region1_id, region2_id, is_final_split)){
        return false;
    }

    return true;

}

bool ScoringFunction::new_split_ok(
    Plan const &plan, 
    int const region1_id, int const region2_id, 
    bool const is_final_split
) const{
    if(!any_hard_constraints) return true;

    auto const num_regions = plan.num_regions;


    // check if new regions are multidistricts 
    bool const is_region1_district = map_params.is_district[plan.region_sizes[region1_id]];
    bool const is_region2_district = map_params.is_district[plan.region_sizes[region2_id]];

    // check if any hard constraints false then immediately return false
    for(auto const &region_constraint_ptr: region_constraint_ptrs){
        // skip if not a hard constraint 
        if(!region_constraint_ptr->hard_constraint) continue;

        // only check if its a district or we're scoring multidistricts 
        if(is_region1_district || !region_constraint_ptr->score_districts_only){
            if(!region_constraint_ptr->region_constraint_ok(plan, region1_id)) return false;
        }
        if(is_region2_district || !region_constraint_ptr->score_districts_only){
            if(!region_constraint_ptr->region_constraint_ok(plan, region2_id)) return false;
        }
    }
    for(auto const &plan_constraint_ptr: plan_constraint_ptrs){
        // skip if not a hard constraint 
        if(!plan_constraint_ptr->hard_constraint) continue;
        // skip if not final plan and only score final plans 
        if(!plan_constraint_ptr->num_regions_to_score[num_regions]) continue;

        if(!plan_constraint_ptr->plan_constraint_ok(plan)) return false;
    }    

    // if final then return 
    if(is_final_split) return true;

    // This is unneccesary right now as the only region constraint not 
    // applied on the final round is pop temper but that doesn't support
    // thresholding right now but this might be useful in the future 
    for(auto const &region_constraint_ptr: non_final_region_constraint_ptrs){
        // skip if not a hard constraint 
        if(!region_constraint_ptr->hard_constraint) continue;

        // only check if its a district or we're scoring multidistricts 
        if(is_region1_district || !region_constraint_ptr->score_districts_only){
            if(!region_constraint_ptr->region_constraint_ok(plan, region1_id)) return false;
        }
        if(is_region2_district || !region_constraint_ptr->score_districts_only){
            if(!region_constraint_ptr->region_constraint_ok(plan, region2_id)) return false;
        }
    }

    return true;
    
}


double ScoringFunction::compute_full_split_plan_soft_score(int const num_regions, 
    PlanVector const &region_ids, RegionSizes const &region_sizes, IntPlanAttribute const &region_pops,
    int const split_region1, int const split_region2
) const{
    double score = 0.0;
    // check if new regions are multidistricts 
    bool const is_region1_district = map_params.is_district[region_sizes[split_region1]];
    bool const is_region2_district = map_params.is_district[region_sizes[split_region2]];


    for(auto const &region_constraint_ptr: region_constraint_ptrs){
        // only check if its a district or we're scoring multidistricts 
        if(is_region1_district || !region_constraint_ptr->score_districts_only){
            score += region_constraint_ptr->compute_soft_region_score(
            num_regions, region_ids, region_sizes, region_pops,
            split_region1);
        }
        if(is_region2_district || !region_constraint_ptr->score_districts_only){
            score += region_constraint_ptr->compute_soft_region_score(
            num_regions, region_ids, region_sizes, region_pops,
            split_region2);
        }
    }
    for(auto const &plan_constraint_ptr: plan_constraint_ptrs){
        // skip if not scoring this many regions 
        if(!plan_constraint_ptr->num_regions_to_score[num_regions]) continue;

        score += plan_constraint_ptr->compute_plan_score(
            num_regions, region_ids, region_sizes, region_pops
        ).second;
    }    


    // This is unneccesary right now as the only region constraint not 
    // applied on the final round is pop temper but that doesn't support
    // thresholding right now but this might be useful in the future 
    for(auto const &region_constraint_ptr: non_final_region_constraint_ptrs){
        // only check if its a district or we're scoring multidistricts 
        if(is_region1_district || !region_constraint_ptr->score_districts_only){
            score += region_constraint_ptr->compute_soft_region_score(
            num_regions, region_ids, region_sizes, region_pops,
            split_region1);
        }
        if(is_region2_district || !region_constraint_ptr->score_districts_only){
            score += region_constraint_ptr->compute_soft_region_score(
            num_regions, region_ids, region_sizes, region_pops,
            split_region2);
        }
    }

    return score;
}