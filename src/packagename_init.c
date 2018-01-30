//
//  packagename_init.c
//  
//
//  Created by Tianying on 1/30/18.
//
//

#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void solveProj_withS(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"solveProj_withS", (DL_FUNC) &solveProj_withS, 10},
    {NULL, NULL, 0}
};

void R_init_DAP(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
