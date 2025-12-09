#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _dcortools_aijbijmem(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dcortools_aijbijmemvec(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dcortools_colsums_subset(SEXP, SEXP);
extern SEXP _dcortools_dcovtermsmem(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dcortools_dcovtermsmemvec(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dcortools_dvartermsmem(SEXP, SEXP, SEXP);
extern SEXP _dcortools_dvartermsmemvec(SEXP, SEXP, SEXP);
extern SEXP _dcortools_hadamard_product(SEXP, SEXP);
extern SEXP _dcortools_matrix_prod_sum(SEXP, SEXP);
extern SEXP _dcortools_matrix_prod_sum_sample(SEXP, SEXP, SEXP);
extern SEXP _dcortools_matrix_prod_sum_subset(SEXP, SEXP, SEXP);
extern SEXP _dcortools_matrix_prod_sum2(SEXP, SEXP);
extern SEXP _dcortools_matrix_sum(SEXP);
extern SEXP _dcortools_normalize_matrix(SEXP, SEXP, SEXP);
extern SEXP _dcortools_specific_vector_prod_sum(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dcortools_squaresym(SEXP);
extern SEXP _dcortools_sum_hadamard_power2(SEXP);
extern SEXP _dcortools_sum_hadamard_power3(SEXP);
extern SEXP _dcortools_SUMAIJBIJ(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dcortools_vector_prod_sum(SEXP, SEXP);
extern SEXP _dcortools_vector_prod_sum_sample(SEXP, SEXP, SEXP);
extern SEXP _dcortools_vector_product(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_dcortools_aijbijmem",                (DL_FUNC) &_dcortools_aijbijmem,                 6},
    {"_dcortools_aijbijmemvec",             (DL_FUNC) &_dcortools_aijbijmemvec,              6},
    {"_dcortools_colsums_subset",           (DL_FUNC) &_dcortools_colsums_subset,            2},
    {"_dcortools_dcovtermsmem",             (DL_FUNC) &_dcortools_dcovtermsmem,              8},
    {"_dcortools_dcovtermsmemvec",          (DL_FUNC) &_dcortools_dcovtermsmemvec,           8},
    {"_dcortools_dvartermsmem",             (DL_FUNC) &_dcortools_dvartermsmem,              3},
    {"_dcortools_dvartermsmemvec",          (DL_FUNC) &_dcortools_dvartermsmemvec,           3},
    {"_dcortools_hadamard_product",         (DL_FUNC) &_dcortools_hadamard_product,          2},
    {"_dcortools_matrix_prod_sum",          (DL_FUNC) &_dcortools_matrix_prod_sum,           2},
    {"_dcortools_matrix_prod_sum_sample",   (DL_FUNC) &_dcortools_matrix_prod_sum_sample,    3},
    {"_dcortools_matrix_prod_sum_subset",   (DL_FUNC) &_dcortools_matrix_prod_sum_subset,    3},
    {"_dcortools_matrix_prod_sum2",         (DL_FUNC) &_dcortools_matrix_prod_sum2,          2},
    {"_dcortools_matrix_sum",               (DL_FUNC) &_dcortools_matrix_sum,                1},
    {"_dcortools_normalize_matrix",         (DL_FUNC) &_dcortools_normalize_matrix,          3},
    {"_dcortools_specific_vector_prod_sum", (DL_FUNC) &_dcortools_specific_vector_prod_sum,  6},
    {"_dcortools_squaresym",                (DL_FUNC) &_dcortools_squaresym,                 1},
    {"_dcortools_sum_hadamard_power2",      (DL_FUNC) &_dcortools_sum_hadamard_power2,       1},
    {"_dcortools_sum_hadamard_power3",      (DL_FUNC) &_dcortools_sum_hadamard_power3,       1},
    {"_dcortools_SUMAIJBIJ",                (DL_FUNC) &_dcortools_SUMAIJBIJ,                10},
    {"_dcortools_vector_prod_sum",          (DL_FUNC) &_dcortools_vector_prod_sum,           2},
    {"_dcortools_vector_prod_sum_sample",   (DL_FUNC) &_dcortools_vector_prod_sum_sample,    3},
    {"_dcortools_vector_product",           (DL_FUNC) &_dcortools_vector_product,            2},
    {NULL, NULL, 0}
};

void R_init_dcortools(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
