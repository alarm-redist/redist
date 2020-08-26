#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _redist_calc_polsbypopper(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_calcPWDh(SEXP);
extern SEXP _redist_countpartitions(SEXP);
extern SEXP _redist_cppGeneratePartitions(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_findBoundary(SEXP, SEXP);
extern SEXP _redist_genAlConn(SEXP, SEXP);
extern SEXP _redist_log_st_map(SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_max_dev(SEXP, SEXP, SEXP);
extern SEXP _redist_n_removed(SEXP, SEXP, SEXP);
extern SEXP _redist_rsg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_sample_partition(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_segregationcalc(SEXP, SEXP, SEXP);
extern SEXP _redist_smc_plans(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _redist_swMH(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_redist_calc_polsbypopper",     (DL_FUNC) &_redist_calc_polsbypopper,      7},
    {"_redist_calcPWDh",              (DL_FUNC) &_redist_calcPWDh,               1},
    {"_redist_countpartitions",       (DL_FUNC) &_redist_countpartitions,        1},
    {"_redist_cppGeneratePartitions", (DL_FUNC) &_redist_cppGeneratePartitions,  7},
    {"_redist_findBoundary",          (DL_FUNC) &_redist_findBoundary,           2},
    {"_redist_genAlConn",             (DL_FUNC) &_redist_genAlConn,              2},
    {"_redist_log_st_map",            (DL_FUNC) &_redist_log_st_map,             4},
    {"_redist_max_dev",               (DL_FUNC) &_redist_max_dev,                3},
    {"_redist_n_removed",             (DL_FUNC) &_redist_n_removed,              3},
    {"_redist_rsg",                   (DL_FUNC) &_redist_rsg,                    6},
    {"_redist_sample_partition",      (DL_FUNC) &_redist_sample_partition,       5},
    {"_redist_segregationcalc",       (DL_FUNC) &_redist_segregationcalc,        3},
    {"_redist_smc_plans",             (DL_FUNC) &_redist_smc_plans,             20},
    {"_redist_swMH",                  (DL_FUNC) &_redist_swMH,                  31},
    {NULL, NULL, 0}
};

void R_init_redist(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
