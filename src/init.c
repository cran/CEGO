#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP lexPermOrder(SEXP);
extern SEXP numericDistanceLevenshtein(SEXP, SEXP);
extern SEXP numericDistanceLongestCommonSubstring(SEXP, SEXP);
extern SEXP permutationDistanceAdjacency(SEXP, SEXP);
extern SEXP permutationDistanceInsert(SEXP, SEXP);
extern SEXP permutationDistanceInterchange(SEXP, SEXP);
extern SEXP permutationDistanceLee(SEXP, SEXP);
extern SEXP permutationDistanceLevenshtein(SEXP, SEXP);
extern SEXP permutationDistanceLongestCommonSubstring(SEXP, SEXP);
extern SEXP permutationDistancePosition(SEXP, SEXP);
extern SEXP permutationDistancePosition2(SEXP, SEXP);
extern SEXP permutationDistanceR(SEXP, SEXP);
extern SEXP permutationDistanceSwap(SEXP, SEXP);
extern SEXP stringDistanceHamming(SEXP, SEXP);
extern SEXP stringDistanceLevenshtein(SEXP, SEXP);
extern SEXP stringDistanceLongestCommonSubstring(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"lexPermOrder",                              (DL_FUNC) &lexPermOrder,                              1},
    {"numericDistanceLevenshtein",                (DL_FUNC) &numericDistanceLevenshtein,                2},
    {"numericDistanceLongestCommonSubstring",     (DL_FUNC) &numericDistanceLongestCommonSubstring,     2},
    {"permutationDistanceAdjacency",              (DL_FUNC) &permutationDistanceAdjacency,              2},
    {"permutationDistanceInsert",                 (DL_FUNC) &permutationDistanceInsert,                 2},
    {"permutationDistanceInterchange",            (DL_FUNC) &permutationDistanceInterchange,            2},
    {"permutationDistanceLee",                    (DL_FUNC) &permutationDistanceLee,                    2},
    {"permutationDistanceLevenshtein",            (DL_FUNC) &permutationDistanceLevenshtein,            2},
    {"permutationDistanceLongestCommonSubstring", (DL_FUNC) &permutationDistanceLongestCommonSubstring, 2},
    {"permutationDistancePosition",               (DL_FUNC) &permutationDistancePosition,               2},
    {"permutationDistancePosition2",              (DL_FUNC) &permutationDistancePosition2,              2},
    {"permutationDistanceR",                      (DL_FUNC) &permutationDistanceR,                      2},
    {"permutationDistanceSwap",                   (DL_FUNC) &permutationDistanceSwap,                   2},
    {"stringDistanceHamming",                     (DL_FUNC) &stringDistanceHamming,                     2},
    {"stringDistanceLevenshtein",                 (DL_FUNC) &stringDistanceLevenshtein,                 2},
    {"stringDistanceLongestCommonSubstring",      (DL_FUNC) &stringDistanceLongestCommonSubstring,      2},
    {NULL, NULL, 0}
};

void attribute_visible R_init_CEGO(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
