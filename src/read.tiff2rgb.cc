/* 
 Authors
 Martin Schlather, schlather@cu.lu

 library for modeling water flux and solute transport:
   makes TIFFReadRGBAImage available in R
    
 Copyright (C)  2003 -- 2006  Martin Schlather, 

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

// #ifdef unix

#include <R.h>
#include <Rdefines.h>
#include <assert.h>
#include "SoPhy.h"
#include "win_linux_aux.h"

SEXP writetiff(SEXP picture, SEXP filename) {
  FILE *file;
  SEXP ret;
  int *p, l2, l, i;
  errno = 0;
  ret=allocVector(INTSXP, 1);
  INTEGER(ret)[0] = 0;
  if ((file = fopen(STRING_VALUE(filename), "w"))==NULL)
    goto ErrorHandling;
  PROTECT(picture);
  p = INTEGER_POINTER(picture);
  l = LENGTH(picture) / 3;
  l2= l * 2;
  for (i=0; i<l; i++) {
    fputc((char) p[i], file);
    fputc((char) p[i+l], file);
    fputc((char) p[i+l2], file);
  }
  fclose(file);
  UNPROTECT(1);
  return ret;

 ErrorHandling:
  INTEGER(ret)[0] = errno;
  PRINTF("open file error occured: %s \n",  strerror(INTEGER(ret)[0]));
  return ret;
}

SEXP readtiff(SEXP filename) {
  unsigned long int Length, Width;
  int x, Error, *r;
  SEXP raster;
  SEXP dim;
  char fname[1000];

  assert(sizeof(x)==4);

  strncpy(fname, CHAR(STRING_ELT(filename, 0)), 1000);
  if ((Error = getTIFFinfo(fname, &Length, &Width)) != 0)
    goto ErrorHandling2;

  PROTECT(dim=allocVector(INTSXP, 3));
  INTEGER_POINTER(dim)[0] = Width;
  INTEGER_POINTER(dim)[1] = Length;
  INTEGER_POINTER(dim)[2] = 4;
  PROTECT(raster=allocArray(INTSXP, dim));
  r = (int*) INTEGER_POINTER(raster);

  strncpy(fname, CHAR(STRING_ELT(filename, 0)), 1000);
  if ((Error = getTIFFimage(fname, Length, Width, r)) != 0)
    goto ErrorHandling;
  UNPROTECT(2);
  return raster;

 ErrorHandling:
  UNPROTECT(2);
 ErrorHandling2:
  PRINTF("error nr %d occured. \n", Error);
  raster=allocVector(INTSXP, 0);
  return raster;
}

//#endif // unix
