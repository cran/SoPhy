
/* 
 Authors
 Martin Schlather, schlather@cu.lu

 library for modeling water flux and solute transport:
    fast cumsum for matrices
    
 Copyright (C)  2003 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

extern "C" {
#include <R.h>
#include <Rdefines.h>
#include <assert.h>
#include "SoPhy.h"
  //#include "auxiliary.h"
}

void anyinside(double *xx, int *nrow, int *ncol, double *min, double *max,
	       int *idx) {
  long int i,t,z;
  for (i=0; i<*nrow; i++) idx[i]=0;
  for (z=t=0; t<*ncol; t++) {
    for (i=0; i<*nrow; i++, z++) {
      if (xx[z]>=*min && xx[z]<=*max) idx[i]=1;
    }
  }
}

SEXP matrixCumsum(SEXP matrix, SEXP margin) {
  // ACHTUNG applyCumsum behaeltdie Orientierung der Matrix bei --
  // nicht wie apply, das die Ergebnisvektoren immer cbind-ed
  SEXP res, dim;
  int d1, d2, i, j, endfor, endfor2;
  double *m, *r;
  PROTECT(dim = NEW_INTEGER(2));
  matrix = AS_NUMERIC(matrix);
  PROTECT(matrix);
  dim = GET_DIM(matrix);
  d1 = INTEGER_POINTER(dim)[0];
  d2 = INTEGER_POINTER(dim)[1];
  PROTECT(res = allocMatrix(REALSXP, d1, d2));
  m = &(NUMERIC_POINTER(matrix)[0]);
  r = &(NUMERIC_POINTER(res)[0]);

  
  if (INTEGER_VALUE(margin)==1) { 
    for (i=0; i<d1; i++) {
      endfor = d2 * d1 + i;
      r[i] = m[i]; 
      for (j=i+d1; j<endfor; j+=d1) {
	r[j] = r[j - d1] + m[j];
      }
    }
  } else if (INTEGER_VALUE(margin)==2) { 
    endfor2 = d2 * d1;
    for (i=0; i<endfor2; i+=d1) {
      endfor = i + d1;
      r[i] = m[i]; 
      for (j=i+1; j<endfor; j++) r[j] = r[j - 1] + m[j];
    }
  } else {
    UNPROTECT(3);
    // printf("%d\n", INTEGER_VALUE(margin));
    error("invalid margin number");
  }
  UNPROTECT(3);
  return res;
}
