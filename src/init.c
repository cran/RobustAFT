#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* How was this made?   AR 05/2018
Within the top level of the package source:
tools::package_native_routine_registration_skeleton('c:/data/R/R-3.4.3/library/RobustAFT',,,FALSE)
Copy all text results to file init.c

Replace void * by int*, float* or double* for each subroutine
To avoid duplicated symbols, subroutine comval and dfcomn2 are renamed comvalz and dfcomn2z
*/

/* .Fortran calls */
extern void F77_NAME(av_tmlnf)(double*, double*, int*, int*, int*, double*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
extern void F77_NAME(av_tmlwf)(double*, double*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
extern void F77_NAME(chia)(int*, float*, float*);
extern void F77_NAME(comvalz)(int*, float*, float*, float*, float*, float*, float*, float*, float*, int*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, int*);
extern void F77_NAME(dfcomn2z)(int*, float*, float*, float*, float*, float*, float*, float*, float*, int*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, int*);
extern void F77_NAME(int44)(float*, float*, float*, float*, float*, float*, int*, int*, int*, float*, int*, int*, int*, int*, int*, float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, float*, float*, float*, float*, float*, float*, float*, int*, float*, float*);
extern void F77_NAME(int51)(float*, float*, int*, float*, int*, int*, float*, int*, int*, int*, int*, float*, float*, float*);
extern void F77_NAME(int59)(float*, float*);
extern void F77_NAME(int60)(float*, float*);
extern void F77_NAME(int61)(float*, float*);
extern void F77_NAME(int62)(float*, float*);
extern void F77_NAME(int92)(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);
extern void F77_NAME(intz21)(float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, float*, float*, float*, int*, float*, float*, int*);
//extern void F77_NAME(nrm2)(float*, int*, int*, int*, float*);
extern void F77_NAME(psia)(int*, float*, float*);
extern void F77_NAME(pspa)(int*, float*, float*);
extern void F77_NAME(qn)(float*, int*, float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int*);
extern void F77_NAME(rhoa)(int*, float*, float*);
extern void F77_NAME(sigama)(float*, float*, float*, float*, float*, float*, int*, float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, float*, int*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, int*, int*);
extern void F77_NAME(srchiamm)(int*, double*, double*, int*, double*);
extern void F77_NAME(srf0w)(double*, double*, int*, double*);
extern void F77_NAME(srintmw)(int*, double*, double*, double*, double*, double*);
extern void F77_NAME(srpsiamm)(int*, double*, double*, int*, double*);
extern void F77_NAME(srpspamm)(int*, double*, double*, int*, double*);
extern void F77_NAME(srrhoamm)(int*, double*, double*, int*, double*);
extern void F77_NAME(zdfvals)(int*, float*);

//   {"nrm2",     (DL_FUNC) &F77_NAME(nrm2),      5},   

static const R_FortranMethodDef FortranEntries[] = {
    {"comvalz",  (DL_FUNC) &F77_NAME(comvalz),  24},
    {"dfcomn2z", (DL_FUNC) &F77_NAME(dfcomn2z), 24},
    {"av_tmlnf", (DL_FUNC) &F77_NAME(av_tmlnf), 20},
    {"av_tmlwf", (DL_FUNC) &F77_NAME(av_tmlwf), 21},
    {"chia",     (DL_FUNC) &F77_NAME(chia),      3},
    {"int44",    (DL_FUNC) &F77_NAME(int44),    35},
    {"int51",    (DL_FUNC) &F77_NAME(int51),    14},
    {"int59",    (DL_FUNC) &F77_NAME(int59),     2},
    {"int60",    (DL_FUNC) &F77_NAME(int60),     2},
    {"int61",    (DL_FUNC) &F77_NAME(int61),     2},
    {"int62",    (DL_FUNC) &F77_NAME(int62),     2},
    {"int92",    (DL_FUNC) &F77_NAME(int92),    18},
    {"intz21",   (DL_FUNC) &F77_NAME(intz21),   31},
    {"psia",     (DL_FUNC) &F77_NAME(psia),      3},
    {"pspa",     (DL_FUNC) &F77_NAME(pspa),      3},
    {"qn",       (DL_FUNC) &F77_NAME(qn),       13},
    {"rhoa",     (DL_FUNC) &F77_NAME(rhoa),      3},
    {"sigama",   (DL_FUNC) &F77_NAME(sigama),   35},
    {"srchiamm", (DL_FUNC) &F77_NAME(srchiamm),  5},
    {"srf0w",    (DL_FUNC) &F77_NAME(srf0w),     4},
    {"srintmw",  (DL_FUNC) &F77_NAME(srintmw),   6},
    {"srpsiamm", (DL_FUNC) &F77_NAME(srpsiamm),  5},
    {"srpspamm", (DL_FUNC) &F77_NAME(srpspamm),  5},
    {"srrhoamm", (DL_FUNC) &F77_NAME(srrhoamm),  5},
    {"zdfvals",  (DL_FUNC) &F77_NAME(zdfvals),   2},
    {NULL, NULL, 0}
};

void R_init_RobustAFT(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}


