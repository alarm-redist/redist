// [[Rcpp::depends(redistmetrics)]]
#include "map_calc.h"
#include <redistmetrics.h>







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
    int const region1_id, int const region2_id
){

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
    int const region1_id, int const region2_id
){
    // If 1 county then just log(1) = 0
    if (n_cty == 1) return 0;

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

