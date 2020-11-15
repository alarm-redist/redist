#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _redist_agg_p2d(SEXP, SEXP, SEXP);
extern SEXP _redist_bias(SEXP, SEXP);
extern SEXP _redist_biasatv(SEXP, SEXP, SEXP);
extern SEXP _redist_calc_polsbypopper(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_calcPWDh(SEXP);
extern SEXP _redist_countpartitions(SEXP);
extern SEXP _redist_cppGeneratePartitions(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_crsg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_declination(SEXP, SEXP, SEXP);
extern SEXP _redist_dist_dist_diff(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_dseats(SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_dseatsDVS(SEXP);
extern SEXP _redist_DVS(SEXP, SEXP);
extern SEXP _redist_effgap(SEXP, SEXP, SEXP);
extern SEXP _redist_effgapEP(SEXP, SEXP, SEXP);
extern SEXP _redist_findBoundary(SEXP, SEXP);
extern SEXP _redist_genAlConn(SEXP, SEXP);
extern SEXP _redist_hamming(SEXP, SEXP);
extern SEXP _redist_log_st_map(SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_lopsidedwins(SEXP, SEXP, SEXP);
extern SEXP _redist_max_dev(SEXP, SEXP, SEXP);
extern SEXP _redist_meanmedian(SEXP);
extern SEXP _redist_minkowski(SEXP, SEXP, SEXP);
extern SEXP _redist_n_removed(SEXP, SEXP, SEXP);
extern SEXP _redist_responsiveness(SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_rsg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_sample_partition(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_segregationcalc(SEXP, SEXP, SEXP);
extern SEXP _redist_smc_plans(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_swMH(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_taugap(SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_var_info_mat(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_redist_agg_p2d",               (DL_FUNC) &_redist_agg_p2d,                3},
    {"_redist_bias",                  (DL_FUNC) &_redist_bias,                   2},
    {"_redist_biasatv",               (DL_FUNC) &_redist_biasatv,                3},
    {"_redist_calc_polsbypopper",     (DL_FUNC) &_redist_calc_polsbypopper,      7},
    {"_redist_calcPWDh",              (DL_FUNC) &_redist_calcPWDh,               1},
    {"_redist_countpartitions",       (DL_FUNC) &_redist_countpartitions,        1},
    {"_redist_cppGeneratePartitions", (DL_FUNC) &_redist_cppGeneratePartitions,  7},
    {"_redist_crsg",                  (DL_FUNC) &_redist_crsg,                   9},
    {"_redist_declination",           (DL_FUNC) &_redist_declination,            3},
    {"_redist_dist_dist_diff",        (DL_FUNC) &_redist_dist_dist_diff,         7},
    {"_redist_dseats",                (DL_FUNC) &_redist_dseats,                 4},
    {"_redist_dseatsDVS",             (DL_FUNC) &_redist_dseatsDVS,              1},
    {"_redist_DVS",                   (DL_FUNC) &_redist_DVS,                    2},
    {"_redist_effgap",                (DL_FUNC) &_redist_effgap,                 3},
    {"_redist_effgapEP",              (DL_FUNC) &_redist_effgapEP,               3},
    {"_redist_findBoundary",          (DL_FUNC) &_redist_findBoundary,           2},
    {"_redist_genAlConn",             (DL_FUNC) &_redist_genAlConn,              2},
    {"_redist_hamming",               (DL_FUNC) &_redist_hamming,                2},
    {"_redist_log_st_map",            (DL_FUNC) &_redist_log_st_map,             4},
    {"_redist_lopsidedwins",          (DL_FUNC) &_redist_lopsidedwins,           3},
    {"_redist_max_dev",               (DL_FUNC) &_redist_max_dev,                3},
    {"_redist_meanmedian",            (DL_FUNC) &_redist_meanmedian,             1},
    {"_redist_minkowski",             (DL_FUNC) &_redist_minkowski,              3},
    {"_redist_n_removed",             (DL_FUNC) &_redist_n_removed,              3},
    {"_redist_responsiveness",        (DL_FUNC) &_redist_responsiveness,         4},
    {"_redist_rsg",                   (DL_FUNC) &_redist_rsg,                    6},
    {"_redist_sample_partition",      (DL_FUNC) &_redist_sample_partition,       5},
    {"_redist_segregationcalc",       (DL_FUNC) &_redist_segregationcalc,        3},
    {"_redist_smc_plans",             (DL_FUNC) &_redist_smc_plans,             21},
    {"_redist_swMH",                  (DL_FUNC) &_redist_swMH,                  34},
    {"_redist_taugap",                (DL_FUNC) &_redist_taugap,                 4},
    {"_redist_var_info_mat",          (DL_FUNC) &_redist_var_info_mat,           3},
    {NULL, NULL, 0}
};

void R_init_redist(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}