
#ifndef SoPhy_H
#define SoPhy 1

#include "basic.h"


EXTERN SEXP GetHorizons(SEXP h, SEXP nr);
EXTERN SEXP readtiff(SEXP filename);
EXTERN SEXP writetiff(SEXP picture, SEXP filename);
EXTERN void preproot(double *mat, int *cond, double *uptake, int *lmat,
		     double *beta, int *pl, int *lpl,  int *np, int* sequ,
		     int *diridx, int *neuidx, int *atmidx,
		     double *dirval, double *neuval, double *atmval,
		     double *root, int *atmadd, int *Error);
EXTERN SEXP matrixCumsum(SEXP, SEXP margin); 
EXTERN void anyinside(double *xx, int *nrow, int *ncol, 
		      double *min, double *max, int *idx);
#endif
