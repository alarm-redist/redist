#include "cw_main.h"

#include "cw_forest_walk.h"
#include "cw_proposal.h"
#include "lct_graph_plan_type.h"
#include "redist_alg_helpers.h"
#include "scoring.h"
#include "splitting_schedule_types.h"
#include "ust_sampler.h"

#include <RcppThread.h>
#include <climits>
#include <memory>
#include <vector>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List cyclewalk_plans(int N, int warmup, int thin, int ndists, int total_seats,
                           Rcpp::IntegerVector const &district_seat_sizes,
                           Rcpp::List const &adj_list, arma::uvec const &counties,
                           arma::uvec const &pop, double target, double lower, double upper,
                           double compactness, Rcpp::IntegerMatrix const &init_plan,
                           Rcpp::IntegerMatrix const &init_seats, Rcpp::List const &control,
                           Rcpp::List const &constraints, Rcpp::List const &edge_weights,
                           int instep, double cycle_walk_frac, int verbosity,
                           bool diagnostic_mode) {
    // Validate columns.
    if (init_plan.ncol() != 1 || init_seats.ncol() != 1) {
        throw Rcpp::exception("init_plan and init_seats must each be one-column matrices.");
    }

    SamplingSpace const sampling_space = SamplingSpace::LCTGraphSpace;

    // MapParams: graph, counties, pop, sizes, MMD config.
    MapParams const map_params(adj_list, counties, pop, ndists, total_seats,
                               as<std::vector<int>>(district_seat_sizes), lower, target, upper,
                            sampling_space);
    int const V = map_params.V;
    bool const mmd_plans = map_params.is_mmd;

    

    // Splitting schedule (cyclewalk doesn't split in its inner loop, but the
    // schedule is required by USTSampler for the initial tree draws and to
    // satisfy PlanEnsemble's signature). PureMSSplittingSchedule is fine.
    auto splitting_schedule_ptr = std::make_unique<PureMSSplittingSchedule>(
        ndists, total_seats, as<std::vector<int>>(district_seat_sizes));

    // ScoringFunction: replaces calc_gibbs_tgt.
    ScoringFunction const scoring_function(map_params, constraints,
                                           /*pop_temper=*/0, /*smc=*/false);

    // Per-chain RNG.
    int const seed = static_cast<int>(Rcpp::sample(INT_MAX, 1)[0]);
    RNGState rng_state(seed);
    std::vector<RNGState> rng_states;
    rng_states.reserve(1);
    rng_states.emplace_back(seed, 6);
    // Legacy: keep GLOBAL_RNG consistent with set.seed() for any code paths
    // that still reach it.
    global_seed_rng(static_cast<int>(Rcpp::sample(INT_MAX, 1)[0]));

    // Output buffers.
    int const nsims = N / thin;
    if (nsims < 1)
        throw Rcpp::exception("N must be at least thin (need at least one recorded sample).");

    Rcpp::IntegerMatrix saved_plans_mat(V, nsims);
    Rcpp::IntegerMatrix saved_plan_sizes(mmd_plans ? ndists : 1, mmd_plans ? nsims : 1);
    Rcpp::IntegerVector mh_decisions(nsims);

    // Diagnostic buffers (size up to the total number of inner steps).
    int const total_steps = (N + warmup) * std::max(instep, 1);
    std::vector<double> diag_accept_prob;
    std::vector<int> diag_cycle_length;
    std::vector<int> diag_n_valid_cuts;
    if (diagnostic_mode) {
        diag_accept_prob.reserve(total_steps);
        diag_cycle_length.reserve(total_steps);
        diag_n_valid_cuts.reserve(total_steps);
    }

    int cycle_walk_accept = 0;
    int cycle_walk_reject = 0;
    int cycle_walk_fail_no_adj = 0;
    int cycle_walk_fail_few_edges = 0;
    int cycle_walk_fail_no_path = 0;
    int cycle_walk_fail_no_cuts = 0;
    int forest_walk_success = 0;
    int forest_walk_fail = 0;

    if (verbosity >= 1) {
        Rcpp::Rcout.imbue(std::locale::classic());
        Rcpp::Rcout << "CYCLEWALK MCMC\n";
        Rcpp::Rcout << "Sampling " << N << " " << V << "-unit maps with " << ndists
                    << " districts and population between " << lower << " and " << upper
                    << ".\n";
    }

    RObject bar = R_NilValue;

    {
        // Single-thread pool (cyclewalk runs one chain at a time; parallel
        // chains are handled in R via foreach).
        RcppThread::ThreadPool pool(0);

        // Load initial plan into a PlanEnsemble of size 1.
        PlanEnsemble plan_ensemble = get_plan_ensemble(
            map_params, *splitting_schedule_ptr, ndists, /*nsims=*/1, sampling_space, init_plan,
            init_seats, rng_states, pool, verbosity);

        auto *plan = static_cast<LCTGraphPlan *>(plan_ensemble.plan_ptr_vec[0].get());

        // Edge weights and LCT setup.
        plan->set_edge_weights(map_params, edge_weights);

        USTSampler ust_sampler(map_params, *splitting_schedule_ptr);
        if (verbosity >= 1) {
            Rcpp::Rcout << "Initializing spanning forest...\n";
        }
        if (!plan->init_lct_from_regions(map_params, ust_sampler, rng_state)) {
            throw Rcpp::exception(
                "Failed to construct initial spanning forest; UST sampler exhausted retries.");
        }

        // Sanity-check initial pop bounds.
        auto pop_check = plan->all_region_pops_valid(map_params);
        if (!pop_check.first) {
            throw Rcpp::exception("Initial plan violates population bounds.");
        }

        if (verbosity >= 3) {
            plan->print_state(verbosity);
        }
        if (verbosity >= 1) {
            bar = cli_progress_bar(N + warmup, cli_config(false));
        }

        // ---- Main loop ----
        int saved_col = 0;
        try {
            for (int i = 1; i <= N + warmup; ++i) {
                bool any_cycle_walk_accepted = false;

                for (int step = 0; step < instep; ++step) {
                    int result;
                    if (rng_state.r_unif() < cycle_walk_frac) {
                        CycleWalkDiagnostics diag;
                        result = cycle_walk(*plan, map_params, scoring_function, compactness,
                                            rng_state, diag);
                        if (diagnostic_mode) {
                            diag_accept_prob.push_back(diag.accept_prob);
                            diag_cycle_length.push_back(diag.cycle_length);
                            diag_n_valid_cuts.push_back(diag.n_valid_cuts);
                        }
                        if (result == 1) {
                            cycle_walk_accept++;
                            any_cycle_walk_accepted = true;
                        } else if (result == 0) {
                            cycle_walk_reject++;
                        } else if (result == -1) {
                            cycle_walk_fail_no_adj++;
                        } else if (result == -2) {
                            cycle_walk_fail_few_edges++;
                        } else if (result == -3) {
                            cycle_walk_fail_no_path++;
                        } else if (result == -4) {
                            cycle_walk_fail_no_cuts++;
                        }
                    } else {
                        result = internal_forest_walk(*plan, map_params, rng_state);
                        if (diagnostic_mode) {
                            diag_accept_prob.push_back(NA_REAL);
                            diag_cycle_length.push_back(0);
                            diag_n_valid_cuts.push_back(0);
                        }
                        if (result == 0)
                            forest_walk_success++;
                        else
                            forest_walk_fail++;
                    }
                }

                // Save post-warmup samples at thinning intervals.
                if (i > warmup && (i - warmup) % thin == 0) {
                    auto col = saved_plans_mat.column(saved_col);
                    for (int v = 0; v < V; ++v) {
                        col[v] = static_cast<int>(plan->region_ids[v]) + 1; // 1-indexed
                    }
                    if (mmd_plans) {
                        auto sz_col = saved_plan_sizes.column(saved_col);
                        for (int d = 0; d < ndists; ++d) {
                            sz_col[d] = static_cast<int>(plan->region_sizes[d]);
                        }
                    }
                    mh_decisions[saved_col] = any_cycle_walk_accepted ? 1 : 0;
                    saved_col++;
                }

                if (verbosity >= 1 && CLI_SHOULD_TICK)
                    cli_progress_set(bar, i - 1);
                if (i % 100 == 0)
                    Rcpp::checkUserInterrupt();
                if (saved_col >= nsims)
                    break;
            }
        } catch (Rcpp::internal::InterruptedException &) {
            if (verbosity >= 1)
                cli_progress_done(bar);
            throw;
        }

        if (verbosity >= 1)
            cli_progress_done(bar);
    }

    if (verbosity >= 2) {
        Rcpp::Rcout << "\n[Forest Walk] Success: " << forest_walk_success
                    << ", Fail: " << forest_walk_fail << "\n";
        Rcpp::Rcout << "[Cycle Walk] Accept: " << cycle_walk_accept
                    << ", Reject: " << cycle_walk_reject << "\n";
        int total_fail = cycle_walk_fail_no_adj + cycle_walk_fail_few_edges +
                         cycle_walk_fail_no_path + cycle_walk_fail_no_cuts;
        if (total_fail > 0) {
            Rcpp::Rcout << "[Cycle Walk Failures] Total: " << total_fail;
            if (cycle_walk_fail_few_edges > 0)
                Rcpp::Rcout << ", <2 edges: " << cycle_walk_fail_few_edges;
            if (cycle_walk_fail_no_cuts > 0)
                Rcpp::Rcout << ", no valid cuts: " << cycle_walk_fail_no_cuts;
            if (cycle_walk_fail_no_path > 0)
                Rcpp::Rcout << ", no path: " << cycle_walk_fail_no_path;
            if (cycle_walk_fail_no_adj > 0)
                Rcpp::Rcout << ", no adj: " << cycle_walk_fail_no_adj;
            Rcpp::Rcout << "\n";
        }
    }

    List diagnostics = List::create(
        Named("accept_prob") =
            diagnostic_mode ? NumericVector(diag_accept_prob.begin(), diag_accept_prob.end())
                            : NumericVector(),
        Named("cycle_length") =
            diagnostic_mode ? IntegerVector(diag_cycle_length.begin(), diag_cycle_length.end())
                            : IntegerVector(),
        Named("n_valid_cuts") =
            diagnostic_mode ? IntegerVector(diag_n_valid_cuts.begin(), diag_n_valid_cuts.end())
                            : IntegerVector(),
        Named("failure_modes") = List::create(Named("no_adj") = cycle_walk_fail_no_adj,
                                              Named("few_edges") = cycle_walk_fail_few_edges,
                                              Named("no_path") = cycle_walk_fail_no_path,
                                              Named("no_cuts") = cycle_walk_fail_no_cuts));

    SEXP plan_sizes_out = mmd_plans ? SEXP(saved_plan_sizes) : R_NilValue;
    List out =
        List::create(Named("plans") = saved_plans_mat, Named("plan_sizes") = plan_sizes_out,
                     Named("mhdecisions") = mh_decisions, Named("diagnostics") = diagnostics,
                     Named("algorithm") = "cyclewalk");
    return out;
}
