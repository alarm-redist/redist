# Package index

## Simulation Algorithm Implementations

Core functions to perform a redistricting simulation

- [`redist_smc()`](http://alarm-redist.org/redist/reference/redist_smc.md)
  : SMC Redistricting Sampler (McCartan and Imai 2023)
- [`redist_mergesplit()`](http://alarm-redist.org/redist/reference/redist_mergesplit.md)
  : Merge-Split/Recombination MCMC Redistricting Sampler (Carter et al.
  2019)
- [`redist_mergesplit_parallel()`](http://alarm-redist.org/redist/reference/redist_mergesplit_parallel.md)
  : Parallel Merge-Split/Recombination MCMC Redistricting Sampler
- [`redist_flip()`](http://alarm-redist.org/redist/reference/redist_flip.md)
  : 'Flip' Markov Chain Monte Carlo Redistricting Simulation (Fifield et
  al. 2020)
- [`redist_shortburst()`](http://alarm-redist.org/redist/reference/redist_shortburst.md)
  : Redistricting Optimization through Short Bursts
- [`redist_constr()`](http://alarm-redist.org/redist/reference/redist_constr.md)
  : Set up constraints for sampling
- [`add_constr_status_quo()`](http://alarm-redist.org/redist/reference/constraints.md)
  [`add_constr_grp_pow()`](http://alarm-redist.org/redist/reference/constraints.md)
  [`add_constr_grp_hinge()`](http://alarm-redist.org/redist/reference/constraints.md)
  [`add_constr_grp_inv_hinge()`](http://alarm-redist.org/redist/reference/constraints.md)
  [`add_constr_compet()`](http://alarm-redist.org/redist/reference/constraints.md)
  [`add_constr_incumbency()`](http://alarm-redist.org/redist/reference/constraints.md)
  [`add_constr_splits()`](http://alarm-redist.org/redist/reference/constraints.md)
  [`add_constr_multisplits()`](http://alarm-redist.org/redist/reference/constraints.md)
  [`add_constr_total_splits()`](http://alarm-redist.org/redist/reference/constraints.md)
  [`add_constr_pop_dev()`](http://alarm-redist.org/redist/reference/constraints.md)
  [`add_constr_segregation()`](http://alarm-redist.org/redist/reference/constraints.md)
  [`add_constr_polsby()`](http://alarm-redist.org/redist/reference/constraints.md)
  [`add_constr_fry_hold()`](http://alarm-redist.org/redist/reference/constraints.md)
  [`add_constr_log_st()`](http://alarm-redist.org/redist/reference/constraints.md)
  [`add_constr_edges_rem()`](http://alarm-redist.org/redist/reference/constraints.md)
  [`add_constr_custom()`](http://alarm-redist.org/redist/reference/constraints.md)
  : Sampling constraints
- [`redist.crsg()`](http://alarm-redist.org/redist/reference/redist.crsg.md)
  : Redistricting via Compact Random Seed and Grow Algorithm
- [`redist.mcmc.mpi()`](http://alarm-redist.org/redist/reference/redist.mcmc.mpi.md)
  : MCMC Redistricting Simulator using MPI
- [`redist.rsg()`](http://alarm-redist.org/redist/reference/redist.rsg.md)
  : Redistricting via Random Seed and Grow Algorithm
- [`redist_flip_anneal()`](http://alarm-redist.org/redist/reference/redist_flip_anneal.md)
  : Flip MCMC Redistricting Simulator using Simulated Annealing

## Analysis Functions

Functions for analyzing simulation outputs

- [`redist_plans()`](http://alarm-redist.org/redist/reference/redist_plans.md)
  : A set of redistricting plans
- [`summary(`*`<redist_plans>`*`)`](http://alarm-redist.org/redist/reference/summary.redist_plans.md)
  : Diagnostic information on sampled plans
- [`add_reference()`](http://alarm-redist.org/redist/reference/add_reference.md)
  : Add a reference plan to a set of plans
- [`subset_sampled()`](http://alarm-redist.org/redist/reference/subset_sampled.md)
  [`subset_ref()`](http://alarm-redist.org/redist/reference/subset_sampled.md)
  : Subset to sampled or reference draws
- [`get_plans_matrix()`](http://alarm-redist.org/redist/reference/get_plans_matrix.md)
  [`as.matrix(`*`<redist_plans>`*`)`](http://alarm-redist.org/redist/reference/get_plans_matrix.md)
  : Extract the matrix of district assignments from a redistricting
  simulation
- [`get_plans_weights()`](http://alarm-redist.org/redist/reference/get_plans_weights.md)
  [`weights(`*`<redist_plans>`*`)`](http://alarm-redist.org/redist/reference/get_plans_weights.md)
  : Extract the sampling weights from a redistricting simulation.
- [`number_by()`](http://alarm-redist.org/redist/reference/number_by.md)
  : Renumber districts to match a quantity of interest
- [`pullback()`](http://alarm-redist.org/redist/reference/pullback.md) :
  Pull back plans to unmerged units
- [`match_numbers()`](http://alarm-redist.org/redist/reference/match_numbers.md)
  : Renumber districts to match an existing plan
- [`tally_var()`](http://alarm-redist.org/redist/reference/tally_var.md)
  : Tally a variable by district
- [`classify_plans()`](http://alarm-redist.org/redist/reference/classify_plans.md)
  : Hierarchically classify a set of redistricting plans
- [`compare_plans()`](http://alarm-redist.org/redist/reference/compare_plans.md)
  : Make a comparison between two sets of plans
- [`is_county_split()`](http://alarm-redist.org/redist/reference/is_county_split.md)
  : Identify which counties are split by a plan
- [`last_plan()`](http://alarm-redist.org/redist/reference/last_plan.md)
  : Extract the last plan from a set of plans
- [`min_move_parity()`](http://alarm-redist.org/redist/reference/min_move_parity.md)
  : Calculates Sparse Population Moves to Minimize Population Deviation
- [`plans_diversity()`](http://alarm-redist.org/redist/reference/plans_diversity.md)
  : Calculate the diversity of a set of plans
- [`plot(`*`<redist_classified>`*`)`](http://alarm-redist.org/redist/reference/plot.redist_classified.md)
  : Plot a plan classification
- [`prec_assignment()`](http://alarm-redist.org/redist/reference/prec_assignment.md)
  : Extract the district assignments for a precinct across all simulated
  plans
- [`prec_cooccurrence()`](http://alarm-redist.org/redist/reference/prec_cooccurrence.md)
  : Compute a matrix of precinct co-occurrences
- [`proj_distr()`](http://alarm-redist.org/redist/reference/proj.md)
  [`proj_avg()`](http://alarm-redist.org/redist/reference/proj.md)
  [`proj_contr()`](http://alarm-redist.org/redist/reference/proj.md) :
  Calculate Projective Distributions, Averages, and Contrasts for a
  Summary Statistic
- [`rbind(`*`<redist_plans>`*`)`](http://alarm-redist.org/redist/reference/rbind.redist_plans.md)
  : Combine multiple sets of redistricting plans
- [`distr_compactness()`](http://alarm-redist.org/redist/reference/redist.compactness.md)
  [`redist.compactness()`](http://alarm-redist.org/redist/reference/redist.compactness.md)
  : Calculate compactness measures for a set of plans
- [`competitiveness()`](http://alarm-redist.org/redist/reference/redist.competitiveness.md)
  [`redist.competitiveness()`](http://alarm-redist.org/redist/reference/redist.competitiveness.md)
  : Compute Competitiveness
- [`plan_distances()`](http://alarm-redist.org/redist/reference/redist.distances.md)
  [`redist.distances()`](http://alarm-redist.org/redist/reference/redist.distances.md)
  : Compute Distance between Partitions
- [`redist.district.splits()`](http://alarm-redist.org/redist/reference/redist.district.splits.md)
  : Counts the Number of Counties within a District
- [`group_frac()`](http://alarm-redist.org/redist/reference/redist.group.percent.md)
  [`redist.group.percent()`](http://alarm-redist.org/redist/reference/redist.group.percent.md)
  : Calculate Group Proportion by District
- [`partisan_metrics()`](http://alarm-redist.org/redist/reference/redist.metrics.md)
  [`redist.metrics()`](http://alarm-redist.org/redist/reference/redist.metrics.md)
  : Calculate gerrymandering metrics for a set of plans
- [`redist.multisplits()`](http://alarm-redist.org/redist/reference/redist.multisplits.md)
  : Counts the Number of Counties Split Between 3 or More Districts
- [`muni_splits()`](http://alarm-redist.org/redist/reference/redist.muni.splits.md)
  [`redist.muni.splits()`](http://alarm-redist.org/redist/reference/redist.muni.splits.md)
  : Counts the Number of Municipalities Split Between Districts
- [`redist.parity()`](http://alarm-redist.org/redist/reference/redist.parity.md)
  [`plan_parity()`](http://alarm-redist.org/redist/reference/redist.parity.md)
  : Calculates Maximum Deviation from Population Parity
- [`segregation_index()`](http://alarm-redist.org/redist/reference/redist.segcalc.md)
  [`redist.segcalc()`](http://alarm-redist.org/redist/reference/redist.segcalc.md)
  : Segregation index calculation for MCMC redistricting.
- [`county_splits()`](http://alarm-redist.org/redist/reference/redist.splits.md)
  [`redist.splits()`](http://alarm-redist.org/redist/reference/redist.splits.md)
  : Count County Splits
- [`redist_ci()`](http://alarm-redist.org/redist/reference/redist_ci.md)
  [`redist_smc_ci()`](http://alarm-redist.org/redist/reference/redist_ci.md)
  [`redist_mcmc_ci()`](http://alarm-redist.org/redist/reference/redist_ci.md)
  : Confidence Intervals for SMC and MCMC Estimates

## Setup Helpers

Functions that help prepare data and select constraints

- [`redist_map()`](http://alarm-redist.org/redist/reference/redist_map.md)
  [`as_redist_map()`](http://alarm-redist.org/redist/reference/redist_map.md)
  :

  Create a `redist_map` object.

- [`get_adj()`](http://alarm-redist.org/redist/reference/get_adj.md)
  [`set_adj()`](http://alarm-redist.org/redist/reference/get_adj.md) :

  Get and set the adjacency graph from a `redist_map` object

- [`get_existing()`](http://alarm-redist.org/redist/reference/get_existing.md)
  :

  Extract the existing district assignment from a `redist_map` object

- [`get_pop_tol()`](http://alarm-redist.org/redist/reference/get_pop_tol.md)
  [`set_pop_tol()`](http://alarm-redist.org/redist/reference/get_pop_tol.md)
  :

  Get and set the population tolerance from a `redist_map` object

- [`get_target()`](http://alarm-redist.org/redist/reference/get_target.md)
  :

  Extract the target district population from a `redist_map` object

- [`is_contiguous()`](http://alarm-redist.org/redist/reference/is_contiguous.md)
  :

  Check that a `redist_map` object is contiguous

- [`merge_by()`](http://alarm-redist.org/redist/reference/merge_by.md) :
  Merge map units

- [`plot(`*`<redist_constr>`*`)`](http://alarm-redist.org/redist/reference/plot.redist_constr.md)
  : Visualize constraints

- [`plot(`*`<redist_map>`*`)`](http://alarm-redist.org/redist/reference/plot.redist_map.md)
  :

  Plot a `redist_map`

- [`redist.coarsen.adjacency()`](http://alarm-redist.org/redist/reference/redist.coarsen.adjacency.md)
  : Coarsen Adjacency List

- [`redist.constraint.helper()`](http://alarm-redist.org/redist/reference/redist.constraint.helper.md)
  : Create Constraints for SMC

- [`redist.county.id()`](http://alarm-redist.org/redist/reference/redist.county.id.md)
  : Create County IDs

- [`redist.county.relabel()`](http://alarm-redist.org/redist/reference/redist.county.relabel.md)
  : Relabel Discontinuous Counties

- [`redist.find.target()`](http://alarm-redist.org/redist/reference/redist.find.target.md)
  : Find Majority Minority Remainder

- [`redist.findparams()`](http://alarm-redist.org/redist/reference/redist.findparams.md)
  :

  Run parameter testing for `redist_flip`

- [`freeze()`](http://alarm-redist.org/redist/reference/redist.freeze.md)
  [`redist.freeze()`](http://alarm-redist.org/redist/reference/redist.freeze.md)
  : Freeze Parts of a Map

- [`make_cores()`](http://alarm-redist.org/redist/reference/redist.identify.cores.md)
  [`redist.identify.cores()`](http://alarm-redist.org/redist/reference/redist.identify.cores.md)
  : Identify Cores of a District (Heuristic)

- [`redist.plot.penalty()`](http://alarm-redist.org/redist/reference/redist.plot.penalty.md)
  : (Deprecated) Visualize Group Power Penalty

- [`redist.reduce.adjacency()`](http://alarm-redist.org/redist/reference/redist.reduce.adjacency.md)
  : Reduce Adjacency List

- [`redist.sink.plan()`](http://alarm-redist.org/redist/reference/redist.sink.plan.md)
  : Sink Plans to 1:ndists

- [`redist.subset()`](http://alarm-redist.org/redist/reference/redist.subset.md)
  : Subset a shp

- [`` `*`( ``*`<redist_scorer>`*`)`](http://alarm-redist.org/redist/reference/scorer-arith.md)
  [`` `+`( ``*`<redist_scorer>`*`)`](http://alarm-redist.org/redist/reference/scorer-arith.md)
  [`` `-`( ``*`<redist_scorer>`*`)`](http://alarm-redist.org/redist/reference/scorer-arith.md)
  : Scoring function arithmetic

- [`combine_scorers()`](http://alarm-redist.org/redist/reference/scorer-combine.md)
  [`cbind(`*`<redist_scorer>`*`)`](http://alarm-redist.org/redist/reference/scorer-combine.md)
  : Combine scoring functions

- [`scorer_group_pct()`](http://alarm-redist.org/redist/reference/scorers.md)
  [`scorer_pop_dev()`](http://alarm-redist.org/redist/reference/scorers.md)
  [`scorer_splits()`](http://alarm-redist.org/redist/reference/scorers.md)
  [`scorer_multisplits()`](http://alarm-redist.org/redist/reference/scorers.md)
  [`scorer_frac_kept()`](http://alarm-redist.org/redist/reference/scorers.md)
  [`scorer_polsby_popper()`](http://alarm-redist.org/redist/reference/scorers.md)
  [`scorer_status_quo()`](http://alarm-redist.org/redist/reference/scorers.md)
  :

  Scoring functions for `redist_shortburst`

## Plotting Tools

Functions for creating plots and maps

- [`plot(`*`<redist_map>`*`)`](http://alarm-redist.org/redist/reference/plot.redist_map.md)
  :

  Plot a `redist_map`

- [`plot(`*`<redist_plans>`*`)`](http://alarm-redist.org/redist/reference/plot.redist_plans.md)
  :

  Summary plots for `\link{redist_plans}`

- [`redist.diagplot()`](http://alarm-redist.org/redist/reference/redist.diagplot.md)
  : Diagnostic plotting functionality for MCMC redistricting.

- [`redist.plot.adj()`](http://alarm-redist.org/redist/reference/redist.plot.adj.md)
  : Creates a Graph Overlay

- [`redist.plot.contr_pfdr()`](http://alarm-redist.org/redist/reference/redist.plot.contr_pfdr.md)
  : Plot a Projective Contrast with positive False Discovery Rate (pFDR)
  Control

- [`redist.plot.cores()`](http://alarm-redist.org/redist/reference/redist.plot.cores.md)
  : Plot Cores

- [`redist.plot.distr_qtys()`](http://alarm-redist.org/redist/reference/redist.plot.distr_qtys.md)
  : Plot quantities by district

- [`redist.plot.hist()`](http://alarm-redist.org/redist/reference/redist.plot.hist.md)
  [`hist(`*`<redist_plans>`*`)`](http://alarm-redist.org/redist/reference/redist.plot.hist.md)
  : Plot a histogram of a summary statistic

- [`redist.plot.majmin()`](http://alarm-redist.org/redist/reference/redist.plot.majmin.md)
  : Majority Minority Plots

- [`redist.plot.map()`](http://alarm-redist.org/redist/reference/redist.plot.map.md)
  : Plot a Map

- [`redist.plot.plans()`](http://alarm-redist.org/redist/reference/redist.plot.plans.md)
  : Plot a district assignment

- [`redist.plot.scatter()`](http://alarm-redist.org/redist/reference/redist.plot.scatter.md)
  : Scatter plot of plan summary statistics

- [`redist.plot.trace()`](http://alarm-redist.org/redist/reference/redist.plot.trace.md)
  : Make a traceplot for a summary statistic

- [`redist.plot.varinfo()`](http://alarm-redist.org/redist/reference/redist.plot.varinfo.md)
  : Static Variation of Information Plot

## Data

Data included to help demonstrate capabilities

- [`EPSG`](http://alarm-redist.org/redist/reference/EPSG.md) : EPSG
  Table
- [`fl25`](http://alarm-redist.org/redist/reference/fl25.md) : Florida
  25 Precinct Shape File
- [`fl250`](http://alarm-redist.org/redist/reference/fl250.md) : Florida
  250 Precinct Shape File
- [`fl25_adj`](http://alarm-redist.org/redist/reference/fl25_adj.md) :
  Florida 25 Precinct File
- [`fl25_enum`](http://alarm-redist.org/redist/reference/fl25_enum.md) :
  All Partitions of 25 Precincts into 3 Congressional Districts (No
  Population Constraint)
- [`fl70`](http://alarm-redist.org/redist/reference/fl70.md) : Florida
  70 Precinct Shape File
- [`iowa`](http://alarm-redist.org/redist/reference/iowa.md) : Iowa
  County File

## Post Processing Helpers

Functions that help setup outputs for easier use

- [`redist.combine.mpi()`](http://alarm-redist.org/redist/reference/redist.combine.mpi.md)
  :

  Combine successive runs of `redist.mcmc.mpi`

- [`redist.ipw()`](http://alarm-redist.org/redist/reference/redist.ipw.md)
  : Inverse probability reweighting for MCMC Redistricting

- [`redist.smc_is_ci()`](http://alarm-redist.org/redist/reference/redist.smc_is_ci.md)
  : (Deprecated) Confidence Intervals for Importance Sampling Estimates

- [`redist.uncoarsen()`](http://alarm-redist.org/redist/reference/redist.uncoarsen.md)
  : Uncoarsen a District Matrix

## Enumeration Tools

Functions for more involved enumeration choices

- [`redist.calc.frontier.size()`](http://alarm-redist.org/redist/reference/redist.calc.frontier.size.md)
  : Calculate Frontier Size
- [`redist.enumpart()`](http://alarm-redist.org/redist/reference/redist.enumpart.md)
  : Enumerate All Parititions (Fifield et al. 2020)
- [`redist.init.enumpart()`](http://alarm-redist.org/redist/reference/redist.init.enumpart.md)
  : Initialize enumpart
- [`redist.prep.enumpart()`](http://alarm-redist.org/redist/reference/redist.prep.enumpart.md)
  : Prepares a run of the enumpart algorithm by ordering edges
- [`redist.read.enumpart()`](http://alarm-redist.org/redist/reference/redist.read.enumpart.md)
  : Read Results from enumpart
- [`redist.run.enumpart()`](http://alarm-redist.org/redist/reference/redist.run.enumpart.md)
  : Runs the enumpart algorithm

## Miscellaneous

Other functions

- [`redist-package`](http://alarm-redist.org/redist/reference/redist-package.md)
  [`redist`](http://alarm-redist.org/redist/reference/redist-package.md)
  : redist: Simulation Methods for Legislative Redistricting

- [`avg_by_prec()`](http://alarm-redist.org/redist/reference/avg_by_prec.md)
  : Average a variable by precinct (Deprecated)

- [`get_mh_acceptance_rate()`](http://alarm-redist.org/redist/reference/get_mh_acceptance_rate.md)
  : Extract the Metropolis Hastings Acceptance Rate

- [`get_sampling_info()`](http://alarm-redist.org/redist/reference/get_sampling_info.md)
  : Extract the sampling information from a redistricting simulation

- [`pl()`](http://alarm-redist.org/redist/reference/pl.md) :

  Access the Current
  [`redist_plans()`](http://alarm-redist.org/redist/reference/redist_plans.md)
  Object

- [`plot(`*`<redist_map>`*`)`](http://alarm-redist.org/redist/reference/plot.redist_map.md)
  :

  Plot a `redist_map`

- [`print(`*`<redist_classified>`*`)`](http://alarm-redist.org/redist/reference/print.redist_classified.md)
  : Print redist_classified objects

- [`print(`*`<redist_constr>`*`)`](http://alarm-redist.org/redist/reference/print.redist_constr.md)
  : Generic to print redist_constr

- [`print(`*`<redist_map>`*`)`](http://alarm-redist.org/redist/reference/print.redist_map.md)
  : Generic to print redist_map

- [`print(`*`<redist_plans>`*`)`](http://alarm-redist.org/redist/reference/print.redist_plans.md)
  :

  Print method for `redist_plans`

- [`redist.adjacency()`](http://alarm-redist.org/redist/reference/redist.adjacency.md)
  : Adjacency List functionality for redist

- [`redist.dist.pop.overlap()`](http://alarm-redist.org/redist/reference/redist.dist.pop.overlap.md)
  : Compare the Population Overlap Across Plans at the District Level

- [`redist.plot.wted.adj()`](http://alarm-redist.org/redist/reference/redist.plot.wted.adj.md)
  : Plot Weighted Border Adjacency

- [`redist.prec.pop.overlap()`](http://alarm-redist.org/redist/reference/redist.prec.pop.overlap.md)
  : Compare the Population Overlap Across Plans at the Precinct Level

- [`redist.random.subgraph()`](http://alarm-redist.org/redist/reference/redist.random.subgraph.md)
  : Return a random subgraph of a shape

- [`redist.reorder()`](http://alarm-redist.org/redist/reference/redist.reorder.md)
  : Reorders district numbers

- [`redist.wted.adj()`](http://alarm-redist.org/redist/reference/redist.wted.adj.md)
  : Create Weighted Adjacency Data

- [`redist_quantile_trunc()`](http://alarm-redist.org/redist/reference/redist_quantile_trunc.md)
  : Helper function to truncate importance weights
