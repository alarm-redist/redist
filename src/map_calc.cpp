// [[Rcpp::depends(redistmetrics)]]
#include "map_calc.h"
#include <redistmetrics.h>



/*
 * Compute the logarithm of the graph theoretic length of the boundary between
 * `region1_id` and `region1_id`
 */
double log_graph_boundary(const Graph &g, const arma::subview_col<arma::uword> &region_ids,
                    int const region1_id, int const region2_id,
                    int const num_counties, arma::uvec counties){
    int V = g.size();
    std::set<int> split_counties;

    if(num_counties > 1 && false){
        // first find the counties that are not wholly contained within either district
        for (int v = 0; v < V; v++)
        {
            int v_region = region_ids(v);
            int v_county = counties(v);
            // ignore if not pairs 
            if(v_region != region1_id && v_region != region2_id) continue;

            for(auto const &u: g[v]){
                int u_region = region_ids(u);
                // ignore if not region 2
                if(u_region == v_region){
                    continue;
                }
                int u_county = counties(u);
                // else if region splits county add it to the set 
                if(v_county == u_county){
                    split_counties.insert(v_county);
                }
            }
        }
    }
    bool const check_county_boundaries = split_counties.size() > 0 && false;

    double count = 0; // number of cuttable edges to create eq-pop districts
    for (int v = 0; v < V; v++) {
        int v_region = region_ids(v);
        if (v_region != region1_id) continue; // Only count if starting vertex in region 1
        for (int u : g[v]) {
            int u_region = region_ids(u);
            // ignore if not the right id
            if (region_ids(u) != region2_id){
                continue;
            }else if( // ignore if we care about counties and this is invalid boundary split
                check_county_boundaries &&
                (counties(v) != counties(u)) &&
                split_counties.count(counties(v)) > 0 &&
                split_counties.count(counties(u)) > 0
            ){ 
                continue;
            }
            // otherwise, boundary with v -> u
            count += 1.0;
        }
    }

    return std::log(count);
}




/*
 * Compute the status quo penalty for district `distr`
 */
double eval_sq_entropy(const subview_col<uword> &districts, const uvec &current,
                       int distr, const uvec &pop, int n_distr, int n_current, int V) {
    double accuml = 0;
    for (int j = 1; j <= n_current; j++) { // 1-indexed districts
        double pop_overlap = 0;
        double pop_total = 0;
        for (int k = 0; k < V; k++) {
            if (current[k] != j) continue;
            pop_total += pop[k];

            if (districts[k] == distr)
                pop_overlap += pop[k];
        }
        double frac = pop_overlap / pop_total;
        if (frac > 0)
            accuml += frac * std::log(frac);
    }

    return -accuml / n_distr / std::log(n_current);
}

/*
 * Compute the new, hinge group penalty for district `distr`
 */
double eval_grp_hinge(const subview_col<uword> &districts, int distr,
                      const vec &tgts_grp, const uvec &grp_pop, const uvec &total_pop) {
    uvec idxs = find(districts == distr);
    double frac = ((double) sum(grp_pop(idxs))) / sum(total_pop(idxs));
    // figure out which to compare it to
    double target;
    double diff = 1;
    int n_tgt = tgts_grp.size();
    for (int i = 0; i < n_tgt; i++) {
        double new_diff = std::fabs(tgts_grp[i] - frac);
        if (new_diff <= diff) {
            diff = new_diff;
            target = tgts_grp[i];
        }
    }

    return std::sqrt(std::max(0.0, target - frac));
}

/*
 * Compute the new, hinge group penalty for district `distr`
 */
double eval_grp_inv_hinge(const subview_col<uword> &districts, int distr,
                      const vec &tgts_grp, const uvec &grp_pop, const uvec &total_pop) {
    uvec idxs = find(districts == distr);
    double frac = ((double) sum(grp_pop(idxs))) / sum(total_pop(idxs));
    // figure out which to compare it to
    double target;
    double diff = 1;
    int n_tgt = tgts_grp.size();
    for (int i = 0; i < n_tgt; i++) {
        double new_diff = std::fabs(tgts_grp[i] - frac);
        if (new_diff <= diff) {
            diff = new_diff;
            target = tgts_grp[i];
        }
    }

    return std::sqrt(std::max(0.0, frac - target));
}

/*
 * Compute the power-based group penalty for district `distr`
 */
double eval_grp_pow(const subview_col<uword> &districts, int distr,
                   const uvec &grp_pop, const uvec &total_pop,
                   double tgt_grp, double tgt_other, double pow) {
    uvec idxs = find(districts == distr);
    double frac = ((double) sum(grp_pop(idxs))) / sum(total_pop(idxs));
    return std::pow(std::fabs(frac - tgt_grp) * std::fabs(frac - tgt_other), pow);
}

/*
 * Compute the incumbent-preserving penalty for district `distr`
 */
double eval_inc(const subview_col<uword> &districts, int distr, const uvec &incumbents) {
    int n_inc = incumbents.size();
    double inc_in_distr = -1.0; // first incumbent doesn't count
    for (int i = 0; i < n_inc; i++) {
        if (districts[incumbents[i] - 1] == distr)
            inc_in_distr++;
    }

    if (inc_in_distr < 0.0) {
        inc_in_distr = 0.0;
    }

    return inc_in_distr;
}


// helper function
// calculates districts which appear in each county (but not zeros)
std::vector<std::set<int>> calc_county_dist(const subview_col<uword> &districts,
                                            const uvec &counties, int n_cty,
                                            bool zero_ok) {
    std::vector<std::set<int>> county_dist(n_cty);
    int V = counties.size();
    for (int i = 0; i < n_cty; i++) {
        county_dist[i] = std::set<int>();
    }
    for (int i = 0; i < V; i++) {
        if (zero_ok || districts[i] > 0) {
            county_dist[counties[i]-1].insert(districts[i]);
        }
    }
    return county_dist;
}

/*
 * Compute the county split penalty for district `distr`
 */
double eval_splits(const subview_col<uword> &districts, int distr,
                   const uvec &counties, int n_cty, bool smc) {
    std::vector<std::set<int>> county_dist = calc_county_dist(districts, counties, n_cty, distr == 0);

    int splits = 0;
    for (int i = 0; i < n_cty; i++) {
        int cty_n_distr = county_dist[i].size();
        // for SMC, just count the split when it crosses the threshold
        // for MCMC there is no sequential nature
        bool cond = smc ? cty_n_distr == 2 : cty_n_distr >= 2;
        if (cond) {
            if (smc) {
                auto search = county_dist[i].find(distr);
                if (search != county_dist[i].end()) {
                    splits++;
                }
            } else {
                splits++;
            }
        }
    }

    return splits;
}

/*
 * Compute the county multisplit penalty for district `distr`
 */
double eval_multisplits(const subview_col<uword> &districts, int distr,
                        const uvec &counties, int n_cty, bool smc) {
    std::vector<std::set<int>> county_dist = calc_county_dist(districts, counties, n_cty, distr == 0);

    double splits = 0;
    for (int i = 0; i < n_cty; i++) {
        int cty_n_distr = county_dist[i].size();
        // for SMC, just count the split when it crosses the threshold
        // for MCMC there is no sequential nature
        bool cond = smc ? cty_n_distr == 3 : cty_n_distr >= 3;
        if (cond) {
            if (smc) {
                auto search = county_dist[i].find(distr);
                if (search != county_dist[i].end()) {
                    splits++;
                }
            } else {
                splits++;
            }
        }
    }

    return splits;
}

/*
 * Compute the total splits penalty for district `distr`
 */
double eval_total_splits(const subview_col<uword> &districts, int distr,
                         const uvec &counties, int n_cty, bool smc) {
    std::vector<std::set<int>> county_dist = calc_county_dist(districts, counties, n_cty, distr == 0);

    double splits = 0;
    for (int i = 0; i < n_cty; i++) {
        int cty_n_distr = county_dist[i].size();
        // no over-counting since every split counts
        if (cty_n_distr > 1) {
            if (smc) {
                auto search = county_dist[i].find(distr);
                if (search != county_dist[i].end()) {
                    splits++;
                }
            } else {
                splits++;
            }
        }
    }

    return splits;
}

/*
 * Compute the Polsby Popper penalty for district `distr`
 */
double eval_polsby(const arma::subview_col<arma::uword> &districts, int distr,
                   const ivec &from,
                   const ivec &to,
                   const vec &area,
                   const vec &perimeter) {
    uvec idxs = find(districts == distr);

    double pi4 = 4.0 * 3.14159265;
    double tot_area = sum(area(idxs));

    double tot_perim = 0.0;

    uvec idx = find(to == distr);
    for (int e = 0; e < idx.size(); e++) {
        if(from(idx(e)) == -1) {
            tot_perim += perimeter(idx(e));
        } else {
            if (districts(from(idx(e))) != distr) {
                tot_perim += perimeter(idx(e));
            }
        }
    }

    double dist_peri2 = pow(tot_perim, 2.0);
    return 1.0 - (pi4 * tot_area / dist_peri2);
}

/*
 * Compute the Fryer-Holden penalty for district `distr`
 */
double eval_fry_hold(const subview_col<uword> &districts, int distr,
                     const uvec &total_pop, mat ssdmat, double denominator = 1.0) {
    uvec idxs = find(districts == distr);
    double ssd = 0.0;

    for (int i = 0; i < idxs.size() - 1; i++) {
        for (int k = i + 1; k < idxs.size(); k++) {
            ssd += (double) ssdmat(idxs(i), idxs(k)) * total_pop(idxs(i)) *
                total_pop(idxs(k));
        }
    }

    return ssd / denominator;
}

/*
 * Compute the population penalty for district `distr`
 */
double eval_pop_dev(const subview_col<uword> &districts, int distr,
                       const uvec &total_pop, double parity) {
    uvec idxs = find(districts == distr);
    double pop = sum(total_pop(idxs));

    return std::pow(pop / parity - 1.0, 2.0);
}


/*
 * Compute the segregation penalty for district `distr`
 */
double eval_segregation(const subview_col<uword> &districts, int distr,
                        const uvec &grp_pop, const uvec &total_pop) {

    int T = sum(total_pop);
    double pAll = (double) sum(grp_pop) / T;
    double denom = (double) 2.0 * T * pAll * (1 - pAll);

    uvec idxs = find(districts == distr);
    double grp = sum(grp_pop(idxs));
    double pop = sum(total_pop(idxs));

    return (double)(pop * std::abs((grp / pop) - pAll) / denom);
}

/*
 * Compute the qps penalty for district `distr`
 */
double eval_qps(const subview_col<uword> &districts, int distr,
                const uvec &total_pop, const uvec &cities, int n_city,
                int nd) {

    vec tally(n_city);
    vec pj(n_city);
    vec j(n_city);
    vec sumpj(n_city);

    uvec idxs_d = find(districts == distr);
    double pop = sum(total_pop(idxs_d));

    for (int i = 0; i < n_city; i++){
        uvec idxs = find(cities == (i + 1));
        idxs = arma::intersect(idxs_d, idxs);
        tally(i) = sum(total_pop(idxs));
        if (tally(i) > 0) {
            j(i) += 1;
        }
    }

    pj = tally / pop;
    sumpj = pj * (1.0 -  pj);
    sumpj = sumpj / (double) nd;

    return sum(sumpj) + log(sum(j));
}

/*
 * Compute the log spanning tree penalty for district `distr`
 */
double eval_log_st(const subview_col<uword> &districts, const Graph g,
                   arma::uvec counties, int ndists) {
    return (double) redistmetrics::log_st_map(g, districts, counties, ndists)[0];
}

/*
 * Compute the edges removed penalty for district `distr`
 */
double eval_er(const subview_col<uword> &districts, const Graph g, int ndists) {
    return (double) redistmetrics::n_removed(g, districts, ndists)[0];
}




/*
 * Compute the cooccurence matrix for a set of precincts indexed by `idxs`,
 * given a collection of plans
 */
mat prec_cooccur(umat m, uvec idxs, int ncores) {
    int v = m.n_rows;
    int n = idxs.n_elem;
    mat out(v, v);

    RcppThread::parallelFor(0, v, [&] (int i) {
        out(i, i) = 1;
        for (int j = 0; j < i; j++) {
            double shared = 0;
            for (int k = 0; k < n; k++) {
                shared += m(i, idxs[k]-1) == m(j, idxs[k]-1);
            }
            shared /= n;
            out(i, j) = shared;
            out(j, i) = shared;
        }
    }, ncores);

    return out;
}

/*
 * Compute the percentage of `group` in each district. Asummes `m` is 1-indexed.
 */
NumericMatrix group_pct(
    IntegerMatrix const &plans_mat, 
    vec const &group_pop, vec const &total_pop, 
    int const n_distr, int const num_threads) {
    int V = plans_mat.nrow();
    int num_plans = plans_mat.ncol();

    NumericMatrix grp_distr(n_distr, num_plans);
    NumericMatrix tot_distr(n_distr, num_plans);

    RcppThread::ThreadPool pool(num_threads > 0 ? num_threads : 0);

    pool.parallelFor(0, num_plans, [&] (unsigned int i) {
        for (int j = 0; j < V; j++) {
            int distr = plans_mat(j, i) - 1;
            grp_distr(distr, i) += group_pop[j];
            tot_distr(distr, i) += total_pop[j];
        }
    });

    pool.wait();


    // divide
    pool.parallelFor(0, num_plans, [&] (unsigned int i) {
        for (int j = 0; j < n_distr; j++) {
            grp_distr(j, i) /= tot_distr(j, i);
        }
    });

    pool.wait();

    return grp_distr;
}

/*
 * Compute the percentage of `group` in each district, and return the `k`-th
 * largest such value. Asummes `m` is 1-indexed.
 */
// [[Rcpp::export]]
NumericVector group_pct_top_k(const IntegerMatrix m, const NumericVector group_pop,
                              const NumericVector total_pop, int k, int n_distr) {
    int v = m.nrow();
    int n = m.ncol();
    NumericVector out(n);

    for (int i = 0; i < n; i++) {
         std::vector<double> grp_distr(n_distr, 0.0);
         std::vector<double> tot_distr(n_distr, 0.0);

        for (int j = 0; j < v; j++) {
            int distr = m(j, i) - 1;
            grp_distr[distr] += group_pop[j];
            tot_distr[distr] += total_pop[j];
        }

        for (int j = 0; j < n_distr; j++) {
            grp_distr[j] /= tot_distr[j];
        }

        std::nth_element(grp_distr.begin(), grp_distr.begin() + k - 1,
                         grp_distr.end(), std::greater<double>());

        out[i] = grp_distr[k - 1];
    }

    return out;
}


/*
 * Tally a variable by district.
 */
// TESTED
// NOTE: Maybe can make parallel version of this? Not sure
NumericMatrix pop_tally(IntegerMatrix const &districts, vec const &pop, int const n_distr,
    int const num_threads) {
    int const num_plans = districts.ncol();
    int const V = districts.nrow();

    NumericMatrix tally(n_distr, num_plans);

    // parallel for loop over each plan 
    RcppThread::parallelFor(0, num_plans, [&] (unsigned int i) {
        for (int j = 0; j < V; j++) {
            int d = districts(j, i) - 1; // districts are 1-indexed
            tally(d, i) = tally(d, i) + pop(j);
        }
    }, num_threads > 0 ? num_threads : 0);

    return tally;
}

// We assume that the population deviations are less than 1
// If they are 1 or greater then you cannot uniquely infer sizes
Rcpp::IntegerMatrix infer_region_sizes(
    Rcpp::IntegerMatrix const &region_pops,
    double const lower, double const upper,
    int const total_seats,
    int const num_threads
){
    // 
    int const num_plans = region_pops.ncol() ;
    int const num_regions = region_pops.nrow();

    bool bounds_issues = false;
    // warn if population bounds aren't tight 
    for (int a_size = 2; a_size <= total_seats; a_size++)
    {
        if(upper * (a_size - 1) >= lower * a_size){
            REprintf("WARNING: Population bounds are not tight for size %d and %d\n",
            a_size-1, a_size);
            Rcpp::warning("Population bounds are not tight, inferring a unique number of seats" 
                " may not be possible.\n");
            bounds_issues = true;
        }
    }
    if(bounds_issues){
        REprintf("WARNING: Population bounds are not tight, inferring a unique number of seats" 
                " may not be possible.\n");
    }

    Rcpp::IntegerMatrix region_sizes(num_regions, num_plans);

    // parallel for loop over each plan 
    RcppThread::parallelFor(0, num_plans, [&] (unsigned int i) {
        // loop over each region 
        for (int j = 0; j < num_regions; j++) {
            // get region pop
            auto region_pop = region_pops(j, i);
            int region_size; bool size_selected = false;
            // find the first instance in which the region is in bounds 
            for (int potential_size = 1; potential_size <= total_seats; potential_size++)
            {
                // see if this size works
                if(lower * potential_size <= region_pop && region_pop <= upper * potential_size){
                    region_size = potential_size; 
                    size_selected = true;
                }
            }
            if(!size_selected){
                REprintf("No valid size could be found for Plan %i\n", i+1);
                throw Rcpp::exception("No valid size could be inferred!\n");    
            }

            region_sizes(j, i) = region_size;
        }
    }, num_threads > 0 ? num_threads : 0);
    
    return region_sizes;
}



/*
 * Create the projective distribution of a variable `x`
 */
// [[Rcpp::export]]
NumericMatrix proj_distr_m(IntegerMatrix districts, const arma::vec x,
                           IntegerVector draw_idx, int n_distr) {
    int n = draw_idx.size();
    int V = districts.nrow();

    NumericMatrix out(V, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < V; j++) {
            int idx = draw_idx[i] - 1;
            out(j, i) = x[n_distr*idx + districts(j, idx) - 1];
        }
    }

    return out;
}

/*
 * Compute the maximum deviation from the equal population constraint.
 */
// TESTED
// Could be parallelized as well
NumericVector max_dev(const IntegerMatrix &districts, const arma::vec &pop, int const n_distr,
                      int const num_threads) {
    int const num_plans = districts.ncol();
    double const target_pop = arma::sum(pop) / n_distr;

    NumericMatrix const dev = pop_tally(districts, pop, n_distr, num_threads) / target_pop - 1.0;
    NumericVector res(num_plans);

    RcppThread::parallelFor(0, num_plans, [&] (unsigned int i) {
        for (int j = 0; j < n_distr; j++) {
            // If deviation at district j bigger then record that
            if (std::fabs(dev(j, i)) > res(i)){
                res(i) = std::fabs(dev(j, i));
            }
        }
    }, num_threads > 0 ? num_threads : 0);

    return res;
}

/*
 * Calculate the deviation for cutting at every edge in a spanning tree.
 */
/*
 * Calculate the deviation for cutting at every edge in a spanning tree.
 */
std::vector<double> tree_dev(Tree &ust, int root, const uvec &pop,
                             double total_pop, double target) {
    int V = pop.size();
    std::vector<int> pop_below(V, 0);
    std::vector<int> parent(V);
    tree_pop(ust, root, pop, pop_below, parent);
    // compile a list of candidate edges to cut
    int idx = 0;
    std::vector<double> devs(V-1);
    for (int i = 0; i < V; i++) {
        if (i == root) continue;
        devs.at(idx) = std::min(std::fabs(pop_below.at(i) - target),
                std::fabs(total_pop - pop_below[i] - target)) / target;
        idx++;
    }

    std::sort(devs.begin(), devs.end());

    return devs;
}


/*
 * Column-wise maximum
 */
// [[Rcpp::export]]
NumericVector colmax(const NumericMatrix x) {
    int nrow = x.nrow();
    int ncol = x.ncol();
    NumericVector out(ncol);
    for (int j = 0; j < ncol; j++) {
        double best = x(0, j);
        for (int i = 1; i < nrow; i++) {
            if (x(i, j) > best) {
                best = x(i, j);
            }
        }
        out[j] = best;
    }

    return out;
}
/*
 * Column-wise minimum
 */
// [[Rcpp::export]]
NumericVector colmin(const NumericMatrix x) {
    int nrow = x.nrow();
    int ncol = x.ncol();
    NumericVector out(ncol);
    for (int j = 0; j < ncol; j++) {
        double best = x(0, j);
        for (int i = 1; i < nrow; i++) {
            if (x(i, j) < best) {
                best = x(i, j);
            }
        }
        out[j] = best;
    }

    return out;
}



// computes log number of spanning trees on region intersect county
// In either a region or a merged region 
double compute_log_region_and_county_spanning_tree(
    Graph const &g, const uvec &counties, int const county,
    PlanVector const &region_ids,
    int const region1_id, int const potential_region2_id
){
    // check if region 2 actually different. 
    // If less than zero set to region1
    int const region2_id = potential_region2_id < 0 ? region1_id : potential_region2_id;
    int const V = g.size();
    // number of precincts in this district
    int K = 0;
    std::vector<int> pos(V); // keep track of positions in subgraph
    int start = 0; // where to start loop below, to save time
    for (int i = 0; i < V; i++) {
        pos[i] = K - 1; // minus one because we're dropping 1st row and column
        // Check if in either region and the county
        if ((region_ids[i] == region1_id || region_ids[i] == region2_id) && 
            counties(i) == county) {
            K++;
            if (K == 2) start = i; // start 2nd vertex
        }
    }
    if (K <= 1) return 0;

    mat adj = zeros<mat>(K-1, K-1); // adjacency matrix (minus 1st row and column)
    for (int i = start; i < V; i++) {
        // ignore if not in either region or the county 
        if ((region_ids[i] != region1_id && region_ids[i] != region2_id) 
            || counties(i) != county) continue;

        int prec = pos.at(i);
        if (prec < 0) continue;
        std::vector<int> nbors = g[i];
        int length = nbors.size();
        int degree = 0; // keep track of index within subgraph
        for (int j = 0; j < length; j++) {
            int nbor = nbors[j];
            if ((region_ids[nbor] != region1_id && region_ids[nbor] != region2_id) 
                || counties(nbor) != county) continue;
            degree++;
            if (pos.at(nbor) < 0) continue;
            adj(prec, pos[nbor]) = -1;
        }
        adj(prec, prec) = degree;
    }

    return arma::log_det_sympd(adj);
}


/*
 * Compute the log number of spanning trees for the contracted (ie county level) graph
 */
// TESTED
double compute_log_county_level_spanning_tree(
    Graph const &g, const uvec &counties, int const n_cty,
    PlanVector const &region_ids,
    int const region1_id, int const potential_region2_id
){
    // If 1 county then just log(1) = 0
    if (n_cty == 1) return 0;

    // check if region 2 actually different. 
    // If less than zero set to region1
    int const region2_id = potential_region2_id < 0 ? region1_id : potential_region2_id;
    int const V = g.size();
    // number of counties in this district
    int K = 0;
    std::vector<int> pos(V); // keep track of positions in subgraph
    std::vector<int> seen(n_cty, -2); // county lookup
    int start = 0;
    for (int i = 0; i < V; i++) {
        if (region_ids[i] != region1_id && region_ids[i] != region2_id) continue;

        if (seen[counties(i)-1] < 0) {
            pos.at(i) = K - 1; // minus one because we're dropping 1st row and column
            seen[counties(i)-1] = K;
            K++;
            if (K == 2) start = i; // start 2nd vertex
        } else {
            pos.at(i) = seen.at(counties(i)-1) - 1;
        }
    }
    if (K <= 1) return 0;

    mat adj = zeros<mat>(K-1, K-1); // adjacency matrix (minus 1st row and column)
    for (int i = start; i < V; i++) {
        if (region_ids[i] != region1_id && region_ids[i] != region2_id) continue;

        int cty = pos[i];
        if (cty < 0) continue; // skip 1st row, col
        std::vector<int> nbors = g[i];
        int length = nbors.size();
        for (int j = 0; j < length; j++) {
            int nbor = nbors.at(j);
            if ((region_ids[nbor] != region1_id && region_ids[nbor] != region2_id) ||
                 pos.at(nbor) == cty) continue;
            adj(cty, cty)++;
            if (pos[nbor] < 0) continue; // ignore 1st row and column
            adj(cty, pos[nbor])--;
        }
    }

    return arma::log_det_sympd(adj);
}

// Given a numeric vector of statistics computed on each district this 
// sorts the statistics within each plan.
// the length district_stats must be a multiple of ndists
Rcpp::NumericVector order_district_stats(
    Rcpp::NumericVector const &district_stats, 
    int const ndists,
    int const num_threads
){
    
    if(district_stats.size() % ndists != 0){
        throw Rcpp::exception("The length of the vector of district statistics must be a multiple of the number of districts\n");
    }else if(ndists <= 1){
        throw Rcpp::exception("Number of districts must be at least 2!\n");
    }

    int num_plans = district_stats.size() / ndists;

    NumericVector ordered_district_stats = clone(district_stats);

    RcppThread::parallelFor(0, num_plans, [&] (unsigned int i) {
        // sort each chunk
        int start_index = i * ndists;
        // we don't subtract 1 since the end index is exclusive!!
        int end_index = start_index + ndists;

        std::sort(
            ordered_district_stats.begin() + start_index, 
            ordered_district_stats.begin() + end_index
        );
    }, num_threads > 0 ? num_threads : 0);

    return ordered_district_stats;
}

/**************************
 * Parallel Versions of redistmetric functions for working with very large plans 
 * Can probably remove in production version
 ****************************/

/*
 * Compute the number of edges removed
 * 
 * Parallel version lifted directly from here
 * https://github.com/alarm-redist/redistmetrics/blob/main/src/kirchhoff.cpp
 */
// TESTED
NumericVector parallel_n_removed(const Graph &g, const IntegerMatrix &districts, int const n_distr, 
                                 int const num_threads) {
    int V = g.size();
    int N = districts.ncol();
    NumericVector n_rem(N);

    RcppThread::parallelFor(0, N, [&] (unsigned int n) {
        double removed = 0.0;
        for (int i = 0; i < V; i++) {
          int dist = districts(i, n);
          std::vector<int> nbors = g[i];
          int length = nbors.size();
          for (int j = 0; j < length; j++) {
            if (districts(nbors[j], n) != dist) removed += 1.0;
          }
        }
        n_rem[n] = removed;
    }, num_threads > 0 ? num_threads : 0);
  
    return n_rem/2;
  }



NumericVector parallel_effgap(
    NumericMatrix const &dcounts, NumericMatrix const &rcounts, 
    int const totvote, int const num_threads){
    NumericVector eg(dcounts.ncol());

    NumericMatrix dwaste(dcounts.nrow(), dcounts.ncol());
    NumericMatrix rwaste(rcounts.nrow(), rcounts.ncol());

    int const num_dcount_cols = dcounts.ncol();
    int const num_dcount_rows = dcounts.nrow();

    RcppThread::parallelFor(0, num_dcount_cols, [&] (unsigned int c) {
        for(int r = 0; r <num_dcount_rows; r++){
            int minwin = floor((dcounts(r,c) + rcounts(r,c))/2.0)+1;
            if(dcounts(r,c) > rcounts(r,c)){
                dwaste(r,c) += (dcounts(r,c) - minwin);
                rwaste(r,c) += rcounts(r,c);
            } else{
                dwaste(r,c) += dcounts(r,c);
                rwaste(r,c) += (rcounts(r,c) - minwin);
            }
        }
    }, num_threads > 0 ? num_threads : 0);


  NumericVector netwaste(dcounts.ncol());
  netwaste = colSums(dwaste) - colSums(rwaste);

  for(int i = 0; i < netwaste.size(); i++){
    eg[i] = netwaste[i]/(double)totvote;
  }

  return eg;
}

/*
 * Parallel version, lifted directly from here
 * https://github.com/alarm-redist/redistmetrics/blob/main/src/kirchhoff.cpp
 */
NumericMatrix parallel_agg_p2d(
    IntegerMatrix const &dm, NumericVector const &vote, 
    int const nd, int const num_threads){
    
    NumericMatrix mat = NumericMatrix(nd, dm.ncol());

    int const num_dm_cols = dm.ncol();
    int const num_dm_rows = dm.nrow();

    RcppThread::parallelFor(0, num_dm_cols, [&] (unsigned int j) {
        for(int i = 0; i < num_dm_rows; i++){
            mat(dm(i,j)-1,j) += vote[i];
        }
    }, num_threads > 0 ? num_threads : 0);

    return mat;
}



/*
 * Parallel version, lifted directly from here
 * https://github.com/alarm-redist/redistmetrics/blob/main/src/kirchhoff.cpp
 */
IntegerVector parallel_dseatsDVS(NumericMatrix const &dvs, int const num_threads){
    IntegerVector dseats = IntegerVector(dvs.ncol());
    int const num_dvs_cols = dvs.ncol();
    int const num_dvs_row = dvs.nrow();

    RcppThread::parallelFor(0, num_dvs_cols, [&] (unsigned int c) {
        for(int r = 0; r < num_dvs_row; r++){
            if(dvs(r,c) > .5){
                dseats[c] += 1;
            }
        }
    }, num_threads > 0 ? num_threads : 0);

    return dseats;
}


/*
 * Parallel version, lifted directly from here
 * https://github.com/alarm-redist/redistmetrics/blob/main/src/kirchhoff.cpp
 */
NumericVector parallel_biasatv(
    NumericMatrix const &dvs, 
    double const v, int const nd, int const num_threads){
    NumericVector dshift = (v) - colMeans(dvs);
    NumericVector rshift = (1-v) - colMeans(dvs);
    NumericMatrix dvs_dshift = clone(dvs);
    NumericMatrix dvs_rshift = clone(dvs);


    int const num_dvs_cols = dvs.ncol();
    int const num_dvs_row = dvs.nrow();

    RcppThread::parallelFor(0, num_dvs_cols, [&] (unsigned int c) {
        for(int r = 0; r < num_dvs_row; r++){
            dvs_dshift(r,c) += dshift(c);
            dvs_rshift(r,c) += rshift(c);
        }
    }, num_threads > 0 ? num_threads : 0);

    NumericVector seat_dshift = (NumericVector)parallel_dseatsDVS(dvs_dshift, num_threads)/(double)nd;
    NumericVector seat_rshift = 1.0 - (NumericVector)parallel_dseatsDVS(dvs_rshift, num_threads)/(double)nd;

    return (seat_rshift - seat_dshift)/2;
}




/*
 * Parallel version, lifted directly from here
 * https://github.com/alarm-redist/redistmetrics/blob/main/src/kirchhoff.cpp
 */
NumericMatrix parallelDVS(
    NumericMatrix const &dcounts, NumericMatrix const &rcounts,
    int const num_threads){

    int const num_dcounts_cols = dcounts.ncol();
    int const num_dcounts_row = dcounts.nrow();

    NumericMatrix mat = NumericMatrix(num_dcounts_row, num_dcounts_cols);

    RcppThread::parallelFor(0, num_dcounts_cols, [&] (unsigned int c) {
        for(int r = 0; r < mat.nrow(); r++){
            mat(r,c) = (double)dcounts(r,c)/(dcounts(r,c)+rcounts(r,c));
        }
    }, num_threads > 0 ? num_threads : 0);

    return mat;
  }




/*
 * Parallel version, lifted directly from here
 * https://github.com/alarm-redist/redistmetrics/blob/main/src/kirchhoff.cpp
 */
Rcpp::IntegerVector parallel_splits(
    const IntegerMatrix &dm, const IntegerVector &community,
    int const nd, int const max_split, int const num_threads,
    bool const skip_last){

    IntegerVector ret(dm.ncol());
    int const nc = sort_unique(community).size();
    
    int const num_dm_cols =  dm.ncol();

    int threads_to_use = num_threads > 1 ? num_threads : 0;
    // by column (aka map)

    RcppThread::ThreadPool pool(num_threads > 0 ? num_threads : 0);

    pool.parallelFor(0, num_dm_cols, [&] (unsigned int c) {
        // DO NOT DO RcppThread::parallelFor when using static variables
        // If you do that then the variables will be kept in a global state from R calls 
        // and this will break if you call this function first with a smaller nc
        static thread_local std::vector<std::vector<bool>> seen(nc, std::vector<bool>(nd, false));
        for (int i = 0; i < nc; i++) {
            std::fill(
                seen[i].begin(),
                seen[i].end(), 
                false
            );
        }


        for(int r = 0; r < dm.nrow(); r++){
            seen[community[r] - 1][dm(r, c) - 1] = true;
        }

        int splits = 0;
        int to = nc;
        if (skip_last) {
            to = nc - 1;
        }

        for (int i = 0; i < to; i++) {
            int tot_split = 0;
            for (int j = 0; j < nd; j++) {
                tot_split += seen[i][j];
                if (tot_split > max_split) {
                    splits++;
                    break;
                }
            }
        }

        ret[c] = splits;
    });

    pool.wait();

    
    return ret;
}




NumericMatrix parallel_polsbypopper(IntegerVector const &from,
                           IntegerVector const &to,
                           NumericVector const &area,
                           NumericVector const &perimeter,
                           IntegerMatrix const &dm,
                           int const nd,
                           int const num_threads) {

  NumericMatrix ret(nd, dm.ncol());

  NumericVector zerovec(nd);
  int ne = from.size();
  double pi4 = 4.0*3.14159265;

  RcppThread::ThreadPool pool(num_threads > 0 ? num_threads : 0);

  pool.parallelFor(0, dm.ncol(), [&] (unsigned int c) {
  //for(int c = 0; c < dm.ncol(); c++){
    // set holders to 0 again 
    static thread_local std::vector<double> dist_area(nd);
    static thread_local std::vector<double> dist_peri(nd);

    std::fill(dist_area.begin(), dist_area.end(), 0.0);
    std::fill(dist_peri.begin(), dist_peri.end(), 0.0);

    // Get a vector of areas ~ just sum
    for(int r = 0; r < dm.nrow(); r++){
      dist_area[dm(r,c) - 1] += area(r);
    }
    // Get a vector of perims ~ sum by id'ing borders
    for(int e = 0; e < ne; e++){
      if(from(e) == -1){
        dist_peri[dm(to(e) - 1, c) - 1] += perimeter(e);
      } else {
        if(dm(from(e) - 1, c) != dm(to(e) - 1, c)){
          dist_peri[dm(to(e) - 1, c) - 1] += perimeter(e);
        }
      }
    }


    for (int d = 0; d < nd; d++) {
        ret(d, c) = pi4 * dist_area[d] / std::pow(dist_peri[d], 2.0);
    }

    });


    pool.wait();

  return ret;
}