#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP redist_calcPWDh(SEXP);
extern SEXP redist_countpartitions(SEXP);
extern SEXP redist_cppGeneratePartitions(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP redist_findBoundary(SEXP, SEXP);
extern SEXP redist_genAlConn(SEXP, SEXP);
extern SEXP redist_rsg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP redist_segregationcalc(SEXP, SEXP, SEXP);
extern SEXP redist_swMH(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"redist_calcPWDh",              (DL_FUNC) &redist_calcPWDh,               1},
    {"redist_countpartitions",       (DL_FUNC) &redist_countpartitions,        1},
    {"redist_cppGeneratePartitions", (DL_FUNC) &redist_cppGeneratePartitions,  7},
    {"redist_findBoundary",          (DL_FUNC) &redist_findBoundary,           2},
    {"redist_genAlConn",             (DL_FUNC) &redist_genAlConn,              2},
    {"redist_rsg",                   (DL_FUNC) &redist_rsg,                    6},
    {"redist_segregationcalc",       (DL_FUNC) &redist_segregationcalc,        3},
    {"redist_swMH",                  (DL_FUNC) &redist_swMH,                  24},
    {NULL, NULL, 0}
};

void R_init_redist(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
