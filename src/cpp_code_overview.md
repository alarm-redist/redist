# C++ Redist Code Overview 
# Author: Philip O'Sullivan

## Description
This document provides an overview of the c++ written or edited by Philip O'Sullivan from 2024-2025. As of writing it is probably about half of the c++ code.


# Important Naming convenctions 

    - `ndists` from the R code still means the number of districts
    - `total_seats` is equivalent to `nseats` in the R code. IE `total_seats` is the total number of seats in a map.


# Classes
Most of the code I wrote tries to use object oriented programming to manage things.
    - `MapParams` Class - This class is largely just a wrapper for all the important map related variables that frequently get passed around to functions. By storing it all in this class it helps limit inputs needed to functions as well as store useful information like the number of counties and help avoid needing to recompute it. 
    - `TreeSplitter` Class - This is the abstract base class used for choosing how to split an edge in a spanning tree. Specifically, given a vector of potential `EdgeCuts` it attempts to select one according to some rule. This abstract class makes it easy to add new forward kernels. 
    - `Plan` Class - This is the abstract base class used for representing a plan. At its core it is a lightweight wrapper for a vector mapping vertex ids to the region they are assigned to. For each of the different sampling spaces there are concrete derived classes. In general all the splitting and calculation code is designed to work on an abstract `Plan` object and the sampling space specifics 

    The derived classes are 
        - `GraphPlan` 
        - `ForestPlan`
        - `LinkingEdgePlan`

    This class has methods that require the following classes as inputs 
        - `MapParams`
        - `USTSampler`
        - `PlanMultigraph`
        - `TreeSplitter`

    The derived classes added later
        - `LCTGraphPlan` (cyclewalk) — `GraphPlan` plus a link-cut-tree spanning forest, cross-region edge index, and per-edge weights. See the cyclewalk section below.


# Porting an external sampler into the modular framework

The cyclewalk MCMC sampler (`src/cw_main.cpp`, `src/cw_proposal.cpp`, `src/cw_forest_walk.cpp`, `src/lct_graph_plan_type.{h,cpp}`) was originally written against an older, self-contained API (its own `LCTPartition` state object, bare `GLOBAL_RNG` access, and a separate constraint path via `mcmc_gibbs.h::calc_gibbs_tgt`). It was retrofitted into the modular `MapParams` / `Plan` / `ScoringFunction` / `USTSampler` / `RNGState` architecture in 2026. The same recipe should work for any future sampler that needs to share infrastructure with `redist_smc` / `redist_mergesplit`. The high-level steps:

1. **Decide on a sampling space and Plan subclass.** Inherit from the existing concrete Plan whose adjacency/boundary semantics match. Cyclewalk picks adjacent district pairs by graph-boundary edge count, so it inherits from `GraphPlan`. The new subclass adds any sampler-specific state (cyclewalk adds an LCT, district roots, cross-edge index, and a sparse edge-weight map) and lazy-initializes that state via an `init_*_from_regions(MapParams, USTSampler, RNGState)` method that the entry point calls after `PlanEnsemble` construction. Constructors are inherited with `using GraphPlan::GraphPlan;` so the Plan/`PlanEnsemble` backing-buffer protocol is preserved.

2. **Register the sampling space and PlanEnsemble dispatch.** Add an enum variant to `SamplingSpace` in `redist_types.h`, extend the parser and `to_str` in `redist_types.cpp`, and add a dispatch case in both `PlanEnsemble` constructors in `redist_alg_helpers.cpp`. For a plan type that only needs `GraphPlan`-style construction (no spanning tree drawn at ensemble time), the dispatch case is a one-line `std::make_unique<NewPlan>(...)` mirroring the `GraphPlan` arm. Treat the corresponding `!use_graph_space` guards (e.g. the `rng_states.size()` check) carefully — extend them to whitelist the new space if appropriate.

3. **Rewrite the entry-point `[[Rcpp::export]]` to match the modular signature shape.** Use `merge_split.cpp::ms_plans` as the template. The signature should take `int N, int warmup, int thin, int ndists, int total_seats, Rcpp::IntegerVector district_seat_sizes, ..., Rcpp::IntegerMatrix init_plan, Rcpp::IntegerMatrix init_seats, Rcpp::List control, Rcpp::List constraints, ...`. In the body, construct `MapParams`, a `SplittingSchedule` (use `PureMSSplittingSchedule` if the sampler doesn't actually split — it's still needed by `USTSampler` and `PlanEnsemble`), `ScoringFunction`, a per-chain `RNGState`, a `USTSampler`, and then load the initial plan via `get_plan_ensemble(..., nsims = 1, ...)`. Cast `plan_ensemble.plan_ptr_vec[0]` to the concrete subclass and call its initializer.

4. **Operate on `Plan&` / `MapParams const&` / `ScoringFunction const&` / `RNGState&` throughout.** Convert proposal/move functions to take these by reference rather than reaching for `GLOBAL_RNG` or scanning ad-hoc state. Use `plan.region_ids` (PlanVector view), `plan.region_pops`, `plan.region_sizes`, `map_params.g`, `map_params.pop`, `map_params.V`, `map_params.lower/target/upper`. For population bounds use **size-aware** comparisons (`map_params.lower * region_sizes[d]`, `map_params.upper * region_sizes[d]`) so MMD support comes for free — `Plan::check_region_pop_valid` already does this if you want to delegate.

5. **Replace `mcmc_gibbs.h::calc_gibbs_tgt` with `ScoringFunction`.** Per-region soft penalties → `scoring_function.compute_region_soft_score(plan, region_id, is_final)`. Plan-level → `scoring_function.compute_plan_score(plan)` (returns `(ok, score)` where the bool flags a hard-constraint failure). For MCMC samplers `is_final` is always `true` (the plan is always fully districted). Sample the MH ratio as `(old_soft - new_soft) + (proposal-ratio terms) + (1 - rho) * (old_log_st - new_log_st)`.

6. **Use `Plan::compute_log_region_spanning_trees(map_params, region_id)` for compactness terms.** Do *not* inline `redistmetrics::log_st_distr` / `log_st_contr` loops over counties. The Plan method already sums `compute_log_region_and_county_spanning_tree_eigen_tri` across the region's intersected counties and adds the contracted-multigraph tau term, using the eigen-based implementation which is faster and more numerically stable than the dense `log_det`. Compute `old_log_st` on the live plan before applying the move, apply the move, then compute `new_log_st`; the ratio is `old - new`.

7. **Snapshot for revert.** Plan attributes are views into a flat ensemble buffer. To support reject/undo, snapshot the affected slices into local `std::vector<RegionID>` / `std::vector<int>` containers via `std::vector<RegionID>(plan.region_ids.begin(), plan.region_ids.end())`. On reject, restore with `std::copy(snap.begin(), snap.end(), plan.region_ids.begin())`. Sampler-specific structures (LCT, cross-edge maps) need their own snapshot/restore — for cyclewalk, the cross-edge map is `std::map`-copyable and the LCT is reverted by undoing the link/cut operations in reverse order.

8. **MMD bookkeeping is mostly free.** `MapParams::is_mmd`, `region_sizes`, and the `init_seats` matrix are already threaded through `PlanEnsemble` and the constraint machinery. The entry point just needs to output a `plan_sizes` matrix of shape `ndists × nsims` when `is_mmd` (or `R_NilValue` otherwise). Per-region pop bounds (step 4) plus the SplittingSchedule's `valid_merge_pair_sizes` and `is_district` checks handle the rest. If the sampler's move preserves seat counts (cyclewalk does), no extra logic is needed; if it changes them (a future MMD-aware merge-split flavor), follow `merging.cpp` for how `region_sizes` are updated.

9. **C++ owns warmup, not R.** New convention: `N` (i.e. `nsims`) is the number of *post-warmup* saved samples; the chain runs `warmup + N * thin` total steps. The R wrapper passes `warmup` through and trims nothing. This matches `redist_mergesplit` / `redist_smc`. Update any R-side post-processing that previously trimmed warmup columns.

10. **Beware compile-time gotchas.** Inline `// foo bar baz` comments that wrap across multiple lines *inside* an `[[Rcpp::export]]` parameter list will break `Rcpp::compileAttributes` — it flattens the declaration to one line and the dangling `//` swallows the next parameter. Keep trailing per-parameter comments on a single line. Forward-declare circularly-referenced types in headers (`tree_splitting.h` forward-declares `ScoringFunction` rather than including `scoring.h`) and pull in the full definition in the corresponding `.cpp`.

The cyclewalk port is the worked example of this recipe. `LCTGraphPlan` (in `src/lct_graph_plan_type.{h,cpp}`) shows the Plan subclass; `cw_main.cpp` shows the entry-point shape; `cw_proposal.cpp` / `cw_forest_walk.cpp` show the move-level adaptation; `R/redist_cyclewalk.R` shows the R-side wiring.