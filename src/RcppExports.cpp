// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "redist_types.h"
#include <RcppArmadillo.h>
#include <RcppThread.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// reduce_adj
List reduce_adj(List adj_list, IntegerVector prec_map, int n_keep);
RcppExport SEXP _redist_reduce_adj(SEXP adj_listSEXP, SEXP prec_mapSEXP, SEXP n_keepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type adj_list(adj_listSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type prec_map(prec_mapSEXP);
    Rcpp::traits::input_parameter< int >::type n_keep(n_keepSEXP);
    rcpp_result_gen = Rcpp::wrap(reduce_adj(adj_list, prec_map, n_keep));
    return rcpp_result_gen;
END_RCPP
}
// collapse_adj
Graph collapse_adj(List graph, const arma::uvec& idxs);
RcppExport SEXP _redist_collapse_adj(SEXP graphSEXP, SEXP idxsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type graph(graphSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idxs(idxsSEXP);
    rcpp_result_gen = Rcpp::wrap(collapse_adj(graph, idxs));
    return rcpp_result_gen;
END_RCPP
}
// coarsen_adjacency
List coarsen_adjacency(List adj, IntegerVector groups);
RcppExport SEXP _redist_coarsen_adjacency(SEXP adjSEXP, SEXP groupsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type adj(adjSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type groups(groupsSEXP);
    rcpp_result_gen = Rcpp::wrap(coarsen_adjacency(adj, groups));
    return rcpp_result_gen;
END_RCPP
}
// get_plan_graph
std::vector<std::set<int>> get_plan_graph(List l, int V, IntegerVector plan, int n_distr);
RcppExport SEXP _redist_get_plan_graph(SEXP lSEXP, SEXP VSEXP, SEXP planSEXP, SEXP n_distrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type V(VSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type plan(planSEXP);
    Rcpp::traits::input_parameter< int >::type n_distr(n_distrSEXP);
    rcpp_result_gen = Rcpp::wrap(get_plan_graph(l, V, plan, n_distr));
    return rcpp_result_gen;
END_RCPP
}
// color_graph
IntegerVector color_graph(List l, IntegerVector plan);
RcppExport SEXP _redist_color_graph(SEXP lSEXP, SEXP planSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type l(lSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type plan(planSEXP);
    rcpp_result_gen = Rcpp::wrap(color_graph(l, plan));
    return rcpp_result_gen;
END_RCPP
}
// polsbypopper
NumericMatrix polsbypopper(IntegerVector from, IntegerVector to, NumericVector area, NumericVector perimeter, IntegerMatrix dm, int nd);
RcppExport SEXP _redist_polsbypopper(SEXP fromSEXP, SEXP toSEXP, SEXP areaSEXP, SEXP perimeterSEXP, SEXP dmSEXP, SEXP ndSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type from(fromSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type to(toSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type area(areaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type perimeter(perimeterSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type dm(dmSEXP);
    Rcpp::traits::input_parameter< int >::type nd(ndSEXP);
    rcpp_result_gen = Rcpp::wrap(polsbypopper(from, to, area, perimeter, dm, nd));
    return rcpp_result_gen;
END_RCPP
}
// genAlConn
List genAlConn(List aList, NumericVector cds);
RcppExport SEXP _redist_genAlConn(SEXP aListSEXP, SEXP cdsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type aList(aListSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cds(cdsSEXP);
    rcpp_result_gen = Rcpp::wrap(genAlConn(aList, cds));
    return rcpp_result_gen;
END_RCPP
}
// findBoundary
NumericVector findBoundary(List fullList, List conList);
RcppExport SEXP _redist_findBoundary(SEXP fullListSEXP, SEXP conListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type fullList(fullListSEXP);
    Rcpp::traits::input_parameter< List >::type conList(conListSEXP);
    rcpp_result_gen = Rcpp::wrap(findBoundary(fullList, conList));
    return rcpp_result_gen;
END_RCPP
}
// contiguity
IntegerVector contiguity(List adj, IntegerVector group);
RcppExport SEXP _redist_contiguity(SEXP adjSEXP, SEXP groupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type adj(adjSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group(groupSEXP);
    rcpp_result_gen = Rcpp::wrap(contiguity(adj, group));
    return rcpp_result_gen;
END_RCPP
}
// cores
List cores(List adj, IntegerVector dm, int k, List cd_within_k);
RcppExport SEXP _redist_cores(SEXP adjSEXP, SEXP dmSEXP, SEXP kSEXP, SEXP cd_within_kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type adj(adjSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dm(dmSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< List >::type cd_within_k(cd_within_kSEXP);
    rcpp_result_gen = Rcpp::wrap(cores(adj, dm, k, cd_within_k));
    return rcpp_result_gen;
END_RCPP
}
// update_conncomp
IntegerVector update_conncomp(IntegerVector dm, IntegerVector kvec, List adj);
RcppExport SEXP _redist_update_conncomp(SEXP dmSEXP, SEXP kvecSEXP, SEXP adjSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type dm(dmSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type kvec(kvecSEXP);
    Rcpp::traits::input_parameter< List >::type adj(adjSEXP);
    rcpp_result_gen = Rcpp::wrap(update_conncomp(dm, kvec, adj));
    return rcpp_result_gen;
END_RCPP
}
// crsg
List crsg(List adj_list, NumericVector population, NumericVector area, NumericVector x_center, NumericVector y_center, int Ndistrict, double target_pop, double thresh, int maxiter);
RcppExport SEXP _redist_crsg(SEXP adj_listSEXP, SEXP populationSEXP, SEXP areaSEXP, SEXP x_centerSEXP, SEXP y_centerSEXP, SEXP NdistrictSEXP, SEXP target_popSEXP, SEXP threshSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type adj_list(adj_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type population(populationSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type area(areaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x_center(x_centerSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y_center(y_centerSEXP);
    Rcpp::traits::input_parameter< int >::type Ndistrict(NdistrictSEXP);
    Rcpp::traits::input_parameter< double >::type target_pop(target_popSEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    rcpp_result_gen = Rcpp::wrap(crsg(adj_list, population, area, x_center, y_center, Ndistrict, target_pop, thresh, maxiter));
    return rcpp_result_gen;
END_RCPP
}
// dist_dist_diff
double dist_dist_diff(int p, int i_dist, int j_dist, NumericVector x_center, NumericVector y_center, NumericVector x, NumericVector y);
RcppExport SEXP _redist_dist_dist_diff(SEXP pSEXP, SEXP i_distSEXP, SEXP j_distSEXP, SEXP x_centerSEXP, SEXP y_centerSEXP, SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type i_dist(i_distSEXP);
    Rcpp::traits::input_parameter< int >::type j_dist(j_distSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x_center(x_centerSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y_center(y_centerSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(dist_dist_diff(p, i_dist, j_dist, x_center, y_center, x, y));
    return rcpp_result_gen;
END_RCPP
}
// log_st_map
NumericVector log_st_map(const Graph& g, const arma::umat& districts, const arma::uvec& counties, int n_distr);
RcppExport SEXP _redist_log_st_map(SEXP gSEXP, SEXP districtsSEXP, SEXP countiesSEXP, SEXP n_distrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Graph& >::type g(gSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type districts(districtsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type counties(countiesSEXP);
    Rcpp::traits::input_parameter< int >::type n_distr(n_distrSEXP);
    rcpp_result_gen = Rcpp::wrap(log_st_map(g, districts, counties, n_distr));
    return rcpp_result_gen;
END_RCPP
}
// n_removed
NumericVector n_removed(const Graph& g, const arma::umat& districts, int n_distr);
RcppExport SEXP _redist_n_removed(SEXP gSEXP, SEXP districtsSEXP, SEXP n_distrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Graph& >::type g(gSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type districts(districtsSEXP);
    Rcpp::traits::input_parameter< int >::type n_distr(n_distrSEXP);
    rcpp_result_gen = Rcpp::wrap(n_removed(g, districts, n_distr));
    return rcpp_result_gen;
END_RCPP
}
// countpartitions
int countpartitions(List aList);
RcppExport SEXP _redist_countpartitions(SEXP aListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type aList(aListSEXP);
    rcpp_result_gen = Rcpp::wrap(countpartitions(aList));
    return rcpp_result_gen;
END_RCPP
}
// calcPWDh
NumericMatrix calcPWDh(NumericMatrix x);
RcppExport SEXP _redist_calcPWDh(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(calcPWDh(x));
    return rcpp_result_gen;
END_RCPP
}
// group_pct_top_k
NumericVector group_pct_top_k(const IntegerMatrix m, const NumericVector group_pop, const NumericVector total_pop, int k, int n_distr);
RcppExport SEXP _redist_group_pct_top_k(SEXP mSEXP, SEXP group_popSEXP, SEXP total_popSEXP, SEXP kSEXP, SEXP n_distrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type m(mSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type group_pop(group_popSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type total_pop(total_popSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type n_distr(n_distrSEXP);
    rcpp_result_gen = Rcpp::wrap(group_pct_top_k(m, group_pop, total_pop, k, n_distr));
    return rcpp_result_gen;
END_RCPP
}
// proj_distr_m
NumericMatrix proj_distr_m(IntegerMatrix districts, const arma::vec x, IntegerVector draw_idx, int n_distr);
RcppExport SEXP _redist_proj_distr_m(SEXP districtsSEXP, SEXP xSEXP, SEXP draw_idxSEXP, SEXP n_distrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type districts(districtsSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type draw_idx(draw_idxSEXP);
    Rcpp::traits::input_parameter< int >::type n_distr(n_distrSEXP);
    rcpp_result_gen = Rcpp::wrap(proj_distr_m(districts, x, draw_idx, n_distr));
    return rcpp_result_gen;
END_RCPP
}
// colmax
NumericVector colmax(const NumericMatrix x);
RcppExport SEXP _redist_colmax(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(colmax(x));
    return rcpp_result_gen;
END_RCPP
}
// colmin
NumericVector colmin(const NumericMatrix x);
RcppExport SEXP _redist_colmin(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(colmin(x));
    return rcpp_result_gen;
END_RCPP
}
// prec_cooccur
arma::mat prec_cooccur(arma::umat m, arma::uvec idxs, int ncores);
RcppExport SEXP _redist_prec_cooccur(SEXP mSEXP, SEXP idxsSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::umat >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type idxs(idxsSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(prec_cooccur(m, idxs, ncores));
    return rcpp_result_gen;
END_RCPP
}
// group_pct
NumericMatrix group_pct(arma::umat m, arma::vec group_pop, arma::vec total_pop, int n_distr);
RcppExport SEXP _redist_group_pct(SEXP mSEXP, SEXP group_popSEXP, SEXP total_popSEXP, SEXP n_distrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::umat >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type group_pop(group_popSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type total_pop(total_popSEXP);
    Rcpp::traits::input_parameter< int >::type n_distr(n_distrSEXP);
    rcpp_result_gen = Rcpp::wrap(group_pct(m, group_pop, total_pop, n_distr));
    return rcpp_result_gen;
END_RCPP
}
// pop_tally
NumericMatrix pop_tally(IntegerMatrix districts, arma::vec pop, int n_distr);
RcppExport SEXP _redist_pop_tally(SEXP districtsSEXP, SEXP popSEXP, SEXP n_distrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type districts(districtsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pop(popSEXP);
    Rcpp::traits::input_parameter< int >::type n_distr(n_distrSEXP);
    rcpp_result_gen = Rcpp::wrap(pop_tally(districts, pop, n_distr));
    return rcpp_result_gen;
END_RCPP
}
// max_dev
NumericVector max_dev(const IntegerMatrix districts, const arma::vec pop, int n_distr);
RcppExport SEXP _redist_max_dev(SEXP districtsSEXP, SEXP popSEXP, SEXP n_distrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type districts(districtsSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type pop(popSEXP);
    Rcpp::traits::input_parameter< int >::type n_distr(n_distrSEXP);
    rcpp_result_gen = Rcpp::wrap(max_dev(districts, pop, n_distr));
    return rcpp_result_gen;
END_RCPP
}
// ms_plans
Rcpp::List ms_plans(int N, List l, const arma::uvec init, const arma::uvec& counties, const arma::uvec& pop, int n_distr, double target, double lower, double upper, double rho, List constraints, List control, int k, int thin, int verbosity);
RcppExport SEXP _redist_ms_plans(SEXP NSEXP, SEXP lSEXP, SEXP initSEXP, SEXP countiesSEXP, SEXP popSEXP, SEXP n_distrSEXP, SEXP targetSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP rhoSEXP, SEXP constraintsSEXP, SEXP controlSEXP, SEXP kSEXP, SEXP thinSEXP, SEXP verbositySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< List >::type l(lSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type init(initSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type counties(countiesSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type pop(popSEXP);
    Rcpp::traits::input_parameter< int >::type n_distr(n_distrSEXP);
    Rcpp::traits::input_parameter< double >::type target(targetSEXP);
    Rcpp::traits::input_parameter< double >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< double >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< List >::type constraints(constraintsSEXP);
    Rcpp::traits::input_parameter< List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type verbosity(verbositySEXP);
    rcpp_result_gen = Rcpp::wrap(ms_plans(N, l, init, counties, pop, n_distr, target, lower, upper, rho, constraints, control, k, thin, verbosity));
    return rcpp_result_gen;
END_RCPP
}
// pareto_dominated
LogicalVector pareto_dominated(arma::mat x);
RcppExport SEXP _redist_pareto_dominated(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(pareto_dominated(x));
    return rcpp_result_gen;
END_RCPP
}
// closest_adj_pop
int closest_adj_pop(IntegerVector adj, int i_dist, NumericVector g_prop);
RcppExport SEXP _redist_closest_adj_pop(SEXP adjSEXP, SEXP i_distSEXP, SEXP g_propSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type adj(adjSEXP);
    Rcpp::traits::input_parameter< int >::type i_dist(i_distSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type g_prop(g_propSEXP);
    rcpp_result_gen = Rcpp::wrap(closest_adj_pop(adj, i_dist, g_prop));
    return rcpp_result_gen;
END_RCPP
}
// rint1
Rcpp::IntegerVector rint1(int n, int max);
RcppExport SEXP _redist_rint1(SEXP nSEXP, SEXP maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type max(maxSEXP);
    rcpp_result_gen = Rcpp::wrap(rint1(n, max));
    return rcpp_result_gen;
END_RCPP
}
// runif1
Rcpp::NumericVector runif1(int n, int max);
RcppExport SEXP _redist_runif1(SEXP nSEXP, SEXP maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type max(maxSEXP);
    rcpp_result_gen = Rcpp::wrap(runif1(n, max));
    return rcpp_result_gen;
END_RCPP
}
// resample_lowvar
arma::ivec resample_lowvar(arma::vec wgts);
RcppExport SEXP _redist_resample_lowvar(SEXP wgtsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type wgts(wgtsSEXP);
    rcpp_result_gen = Rcpp::wrap(resample_lowvar(wgts));
    return rcpp_result_gen;
END_RCPP
}
// plan_joint
NumericMatrix plan_joint(IntegerVector m1, IntegerVector m2, NumericVector pop);
RcppExport SEXP _redist_plan_joint(SEXP m1SEXP, SEXP m2SEXP, SEXP popSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type m2(m2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pop(popSEXP);
    rcpp_result_gen = Rcpp::wrap(plan_joint(m1, m2, pop));
    return rcpp_result_gen;
END_RCPP
}
// renumber_matrix
IntegerMatrix renumber_matrix(IntegerMatrix plans, IntegerVector renumb);
RcppExport SEXP _redist_renumber_matrix(SEXP plansSEXP, SEXP renumbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type plans(plansSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type renumb(renumbSEXP);
    rcpp_result_gen = Rcpp::wrap(renumber_matrix(plans, renumb));
    return rcpp_result_gen;
END_RCPP
}
// solve_hungarian
IntegerMatrix solve_hungarian(NumericMatrix costMatrix);
RcppExport SEXP _redist_solve_hungarian(SEXP costMatrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type costMatrix(costMatrixSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_hungarian(costMatrix));
    return rcpp_result_gen;
END_RCPP
}
// rsg
List rsg(List adj_list, NumericVector population, int Ndistrict, double target_pop, double thresh, int maxiter);
RcppExport SEXP _redist_rsg(SEXP adj_listSEXP, SEXP populationSEXP, SEXP NdistrictSEXP, SEXP target_popSEXP, SEXP threshSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type adj_list(adj_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type population(populationSEXP);
    Rcpp::traits::input_parameter< int >::type Ndistrict(NdistrictSEXP);
    Rcpp::traits::input_parameter< double >::type target_pop(target_popSEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    rcpp_result_gen = Rcpp::wrap(rsg(adj_list, population, Ndistrict, target_pop, thresh, maxiter));
    return rcpp_result_gen;
END_RCPP
}
// k_smallest
NumericVector k_smallest(NumericMatrix x, int k);
RcppExport SEXP _redist_k_smallest(SEXP xSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(k_smallest(x, k));
    return rcpp_result_gen;
END_RCPP
}
// k_biggest
NumericVector k_biggest(NumericMatrix x, int k);
RcppExport SEXP _redist_k_biggest(SEXP xSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(k_biggest(x, k));
    return rcpp_result_gen;
END_RCPP
}
// smc_plans
List smc_plans(int N, List l, const arma::uvec& counties, const arma::uvec& pop, int n_distr, double target, double lower, double upper, double rho, arma::umat districts, int n_drawn, int n_steps, List constraints, List control, int verbosity);
RcppExport SEXP _redist_smc_plans(SEXP NSEXP, SEXP lSEXP, SEXP countiesSEXP, SEXP popSEXP, SEXP n_distrSEXP, SEXP targetSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP rhoSEXP, SEXP districtsSEXP, SEXP n_drawnSEXP, SEXP n_stepsSEXP, SEXP constraintsSEXP, SEXP controlSEXP, SEXP verbositySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< List >::type l(lSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type counties(countiesSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type pop(popSEXP);
    Rcpp::traits::input_parameter< int >::type n_distr(n_distrSEXP);
    Rcpp::traits::input_parameter< double >::type target(targetSEXP);
    Rcpp::traits::input_parameter< double >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< double >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type districts(districtsSEXP);
    Rcpp::traits::input_parameter< int >::type n_drawn(n_drawnSEXP);
    Rcpp::traits::input_parameter< int >::type n_steps(n_stepsSEXP);
    Rcpp::traits::input_parameter< List >::type constraints(constraintsSEXP);
    Rcpp::traits::input_parameter< List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< int >::type verbosity(verbositySEXP);
    rcpp_result_gen = Rcpp::wrap(smc_plans(N, l, counties, pop, n_distr, target, lower, upper, rho, districts, n_drawn, n_steps, constraints, control, verbosity));
    return rcpp_result_gen;
END_RCPP
}
// splits
IntegerVector splits(IntegerMatrix dm, IntegerVector community, int nd, int max_split);
RcppExport SEXP _redist_splits(SEXP dmSEXP, SEXP communitySEXP, SEXP ndSEXP, SEXP max_splitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type dm(dmSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type community(communitySEXP);
    Rcpp::traits::input_parameter< int >::type nd(ndSEXP);
    Rcpp::traits::input_parameter< int >::type max_split(max_splitSEXP);
    rcpp_result_gen = Rcpp::wrap(splits(dm, community, nd, max_split));
    return rcpp_result_gen;
END_RCPP
}
// dist_cty_splits
IntegerMatrix dist_cty_splits(IntegerMatrix dm, IntegerVector community, int nd);
RcppExport SEXP _redist_dist_cty_splits(SEXP dmSEXP, SEXP communitySEXP, SEXP ndSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type dm(dmSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type community(communitySEXP);
    Rcpp::traits::input_parameter< int >::type nd(ndSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_cty_splits(dm, community, nd));
    return rcpp_result_gen;
END_RCPP
}
// swMH
List swMH(List aList, NumericVector cdvec, NumericVector popvec, int nsims, List constraints, double eprob, double pct_dist_parity, NumericVector beta_sequence, NumericVector beta_weights, int lambda, double beta, std::string adapt_beta, int adjswap, int exact_mh, int adapt_eprob, int adapt_lambda, int num_hot_steps, int num_annealing_steps, int num_cold_steps, bool verbose);
RcppExport SEXP _redist_swMH(SEXP aListSEXP, SEXP cdvecSEXP, SEXP popvecSEXP, SEXP nsimsSEXP, SEXP constraintsSEXP, SEXP eprobSEXP, SEXP pct_dist_paritySEXP, SEXP beta_sequenceSEXP, SEXP beta_weightsSEXP, SEXP lambdaSEXP, SEXP betaSEXP, SEXP adapt_betaSEXP, SEXP adjswapSEXP, SEXP exact_mhSEXP, SEXP adapt_eprobSEXP, SEXP adapt_lambdaSEXP, SEXP num_hot_stepsSEXP, SEXP num_annealing_stepsSEXP, SEXP num_cold_stepsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type aList(aListSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cdvec(cdvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type popvec(popvecSEXP);
    Rcpp::traits::input_parameter< int >::type nsims(nsimsSEXP);
    Rcpp::traits::input_parameter< List >::type constraints(constraintsSEXP);
    Rcpp::traits::input_parameter< double >::type eprob(eprobSEXP);
    Rcpp::traits::input_parameter< double >::type pct_dist_parity(pct_dist_paritySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_sequence(beta_sequenceSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_weights(beta_weightsSEXP);
    Rcpp::traits::input_parameter< int >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< std::string >::type adapt_beta(adapt_betaSEXP);
    Rcpp::traits::input_parameter< int >::type adjswap(adjswapSEXP);
    Rcpp::traits::input_parameter< int >::type exact_mh(exact_mhSEXP);
    Rcpp::traits::input_parameter< int >::type adapt_eprob(adapt_eprobSEXP);
    Rcpp::traits::input_parameter< int >::type adapt_lambda(adapt_lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type num_hot_steps(num_hot_stepsSEXP);
    Rcpp::traits::input_parameter< int >::type num_annealing_steps(num_annealing_stepsSEXP);
    Rcpp::traits::input_parameter< int >::type num_cold_steps(num_cold_stepsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(swMH(aList, cdvec, popvec, nsims, constraints, eprob, pct_dist_parity, beta_sequence, beta_weights, lambda, beta, adapt_beta, adjswap, exact_mh, adapt_eprob, adapt_lambda, num_hot_steps, num_annealing_steps, num_cold_steps, verbose));
    return rcpp_result_gen;
END_RCPP
}
// tree_pop
int tree_pop(Tree& ust, int vtx, const arma::uvec& pop, std::vector<int>& pop_below, std::vector<int>& parent);
RcppExport SEXP _redist_tree_pop(SEXP ustSEXP, SEXP vtxSEXP, SEXP popSEXP, SEXP pop_belowSEXP, SEXP parentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Tree& >::type ust(ustSEXP);
    Rcpp::traits::input_parameter< int >::type vtx(vtxSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type pop(popSEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type pop_below(pop_belowSEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type parent(parentSEXP);
    rcpp_result_gen = Rcpp::wrap(tree_pop(ust, vtx, pop, pop_below, parent));
    return rcpp_result_gen;
END_RCPP
}
// var_info_vec
NumericVector var_info_vec(IntegerMatrix m, IntegerVector ref, NumericVector pop);
RcppExport SEXP _redist_var_info_vec(SEXP mSEXP, SEXP refSEXP, SEXP popSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type m(mSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ref(refSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pop(popSEXP);
    rcpp_result_gen = Rcpp::wrap(var_info_vec(m, ref, pop));
    return rcpp_result_gen;
END_RCPP
}
// sample_ust
Tree sample_ust(List l, const arma::uvec& pop, double lower, double upper, const arma::uvec& counties, const std::vector<bool> ignore);
RcppExport SEXP _redist_sample_ust(SEXP lSEXP, SEXP popSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP countiesSEXP, SEXP ignoreSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type l(lSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type pop(popSEXP);
    Rcpp::traits::input_parameter< double >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< double >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type counties(countiesSEXP);
    Rcpp::traits::input_parameter< const std::vector<bool> >::type ignore(ignoreSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_ust(l, pop, lower, upper, counties, ignore));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_redist_reduce_adj", (DL_FUNC) &_redist_reduce_adj, 3},
    {"_redist_collapse_adj", (DL_FUNC) &_redist_collapse_adj, 2},
    {"_redist_coarsen_adjacency", (DL_FUNC) &_redist_coarsen_adjacency, 2},
    {"_redist_get_plan_graph", (DL_FUNC) &_redist_get_plan_graph, 4},
    {"_redist_color_graph", (DL_FUNC) &_redist_color_graph, 2},
    {"_redist_polsbypopper", (DL_FUNC) &_redist_polsbypopper, 6},
    {"_redist_genAlConn", (DL_FUNC) &_redist_genAlConn, 2},
    {"_redist_findBoundary", (DL_FUNC) &_redist_findBoundary, 2},
    {"_redist_contiguity", (DL_FUNC) &_redist_contiguity, 2},
    {"_redist_cores", (DL_FUNC) &_redist_cores, 4},
    {"_redist_update_conncomp", (DL_FUNC) &_redist_update_conncomp, 3},
    {"_redist_crsg", (DL_FUNC) &_redist_crsg, 9},
    {"_redist_dist_dist_diff", (DL_FUNC) &_redist_dist_dist_diff, 7},
    {"_redist_log_st_map", (DL_FUNC) &_redist_log_st_map, 4},
    {"_redist_n_removed", (DL_FUNC) &_redist_n_removed, 3},
    {"_redist_countpartitions", (DL_FUNC) &_redist_countpartitions, 1},
    {"_redist_calcPWDh", (DL_FUNC) &_redist_calcPWDh, 1},
    {"_redist_group_pct_top_k", (DL_FUNC) &_redist_group_pct_top_k, 5},
    {"_redist_proj_distr_m", (DL_FUNC) &_redist_proj_distr_m, 4},
    {"_redist_colmax", (DL_FUNC) &_redist_colmax, 1},
    {"_redist_colmin", (DL_FUNC) &_redist_colmin, 1},
    {"_redist_prec_cooccur", (DL_FUNC) &_redist_prec_cooccur, 3},
    {"_redist_group_pct", (DL_FUNC) &_redist_group_pct, 4},
    {"_redist_pop_tally", (DL_FUNC) &_redist_pop_tally, 3},
    {"_redist_max_dev", (DL_FUNC) &_redist_max_dev, 3},
    {"_redist_ms_plans", (DL_FUNC) &_redist_ms_plans, 15},
    {"_redist_pareto_dominated", (DL_FUNC) &_redist_pareto_dominated, 1},
    {"_redist_closest_adj_pop", (DL_FUNC) &_redist_closest_adj_pop, 3},
    {"_redist_rint1", (DL_FUNC) &_redist_rint1, 2},
    {"_redist_runif1", (DL_FUNC) &_redist_runif1, 2},
    {"_redist_resample_lowvar", (DL_FUNC) &_redist_resample_lowvar, 1},
    {"_redist_plan_joint", (DL_FUNC) &_redist_plan_joint, 3},
    {"_redist_renumber_matrix", (DL_FUNC) &_redist_renumber_matrix, 2},
    {"_redist_solve_hungarian", (DL_FUNC) &_redist_solve_hungarian, 1},
    {"_redist_rsg", (DL_FUNC) &_redist_rsg, 6},
    {"_redist_k_smallest", (DL_FUNC) &_redist_k_smallest, 2},
    {"_redist_k_biggest", (DL_FUNC) &_redist_k_biggest, 2},
    {"_redist_smc_plans", (DL_FUNC) &_redist_smc_plans, 15},
    {"_redist_splits", (DL_FUNC) &_redist_splits, 4},
    {"_redist_dist_cty_splits", (DL_FUNC) &_redist_dist_cty_splits, 3},
    {"_redist_swMH", (DL_FUNC) &_redist_swMH, 20},
    {"_redist_tree_pop", (DL_FUNC) &_redist_tree_pop, 5},
    {"_redist_var_info_vec", (DL_FUNC) &_redist_var_info_vec, 3},
    {"_redist_sample_ust", (DL_FUNC) &_redist_sample_ust, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_redist(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
