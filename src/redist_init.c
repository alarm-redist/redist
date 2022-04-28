#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _redist_calcPWDh(SEXP);
extern SEXP _redist_closest_adj_pop(SEXP, SEXP, SEXP);
extern SEXP _redist_coarsen_adjacency(SEXP, SEXP);
extern SEXP _redist_collapse_adj(SEXP, SEXP);
extern SEXP _redist_colmax(SEXP);
extern SEXP _redist_colmin(SEXP);
extern SEXP _redist_color_graph(SEXP, SEXP);
extern SEXP _redist_contiguity(SEXP, SEXP);
extern SEXP _redist_cores(SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_countpartitions(SEXP);
extern SEXP _redist_crsg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_dist_cty_splits(SEXP, SEXP, SEXP);
extern SEXP _redist_dist_dist_diff(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_findBoundary(SEXP, SEXP);
extern SEXP _redist_genAlConn(SEXP, SEXP);
extern SEXP _redist_get_plan_graph(SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_group_pct(SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_group_pct_top_k(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_k_biggest(SEXP, SEXP);
extern SEXP _redist_k_smallest(SEXP, SEXP);
extern SEXP _redist_log_st_map(SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_max_dev(SEXP, SEXP, SEXP);
extern SEXP _redist_ms_plans(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_n_removed(SEXP, SEXP, SEXP);
extern SEXP _redist_plan_joint(SEXP, SEXP, SEXP);
extern SEXP _redist_polsbypopper(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_pop_tally(SEXP, SEXP, SEXP);
extern SEXP _redist_prec_cooccur(SEXP, SEXP, SEXP);
extern SEXP _redist_reduce_adj(SEXP, SEXP, SEXP);
extern SEXP _redist_reindex(SEXP, SEXP);
extern SEXP _redist_renumber_matrix(SEXP, SEXP);
extern SEXP _redist_rsg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_smc_plans(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_solve_hungarian(SEXP);
extern SEXP _redist_splits(SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_swMH(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_update_conncomp(SEXP, SEXP, SEXP);
extern SEXP _redist_var_info_vec(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_redist_calcPWDh",          (DL_FUNC) &_redist_calcPWDh,           1},
    {"_redist_closest_adj_pop",   (DL_FUNC) &_redist_closest_adj_pop,    3},
    {"_redist_coarsen_adjacency", (DL_FUNC) &_redist_coarsen_adjacency,  2},
    {"_redist_collapse_adj",      (DL_FUNC) &_redist_collapse_adj,       2},
    {"_redist_colmax",            (DL_FUNC) &_redist_colmax,             1},
    {"_redist_colmin",            (DL_FUNC) &_redist_colmin,             1},
    {"_redist_color_graph",       (DL_FUNC) &_redist_color_graph,        2},
    {"_redist_contiguity",        (DL_FUNC) &_redist_contiguity,         2},
    {"_redist_cores",             (DL_FUNC) &_redist_cores,              4},
    {"_redist_countpartitions",   (DL_FUNC) &_redist_countpartitions,    1},
    {"_redist_crsg",              (DL_FUNC) &_redist_crsg,               9},
    {"_redist_dist_cty_splits",   (DL_FUNC) &_redist_dist_cty_splits,    3},
    {"_redist_dist_dist_diff",    (DL_FUNC) &_redist_dist_dist_diff,     7},
    {"_redist_findBoundary",      (DL_FUNC) &_redist_findBoundary,       2},
    {"_redist_genAlConn",         (DL_FUNC) &_redist_genAlConn,          2},
    {"_redist_get_plan_graph",    (DL_FUNC) &_redist_get_plan_graph,     4},
    {"_redist_group_pct",         (DL_FUNC) &_redist_group_pct,          4},
    {"_redist_group_pct_top_k",   (DL_FUNC) &_redist_group_pct_top_k,    5},
    {"_redist_k_biggest",         (DL_FUNC) &_redist_k_biggest,          2},
    {"_redist_k_smallest",        (DL_FUNC) &_redist_k_smallest,         2},
    {"_redist_log_st_map",        (DL_FUNC) &_redist_log_st_map,         4},
    {"_redist_max_dev",           (DL_FUNC) &_redist_max_dev,            3},
    {"_redist_ms_plans",          (DL_FUNC) &_redist_ms_plans,          15},
    {"_redist_n_removed",         (DL_FUNC) &_redist_n_removed,          3},
    {"_redist_plan_joint",        (DL_FUNC) &_redist_plan_joint,         3},
    {"_redist_polsbypopper",      (DL_FUNC) &_redist_polsbypopper,       6},
    {"_redist_pop_tally",         (DL_FUNC) &_redist_pop_tally,          3},
    {"_redist_prec_cooccur",      (DL_FUNC) &_redist_prec_cooccur,       3},
    {"_redist_reduce_adj",        (DL_FUNC) &_redist_reduce_adj,         3},
    {"_redist_reindex",           (DL_FUNC) &_redist_reindex,            2},
    {"_redist_renumber_matrix",   (DL_FUNC) &_redist_renumber_matrix,    2},
    {"_redist_rsg",               (DL_FUNC) &_redist_rsg,                6},
    {"_redist_smc_plans",         (DL_FUNC) &_redist_smc_plans,         15},
    {"_redist_solve_hungarian",   (DL_FUNC) &_redist_solve_hungarian,    1},
    {"_redist_splits",            (DL_FUNC) &_redist_splits,             4},
    {"_redist_swMH",              (DL_FUNC) &_redist_swMH,              20},
    {"_redist_update_conncomp",   (DL_FUNC) &_redist_update_conncomp,    3},
    {"_redist_var_info_vec",      (DL_FUNC) &_redist_var_info_vec,       3},
    {NULL, NULL, 0}
};

void R_init_redist(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
