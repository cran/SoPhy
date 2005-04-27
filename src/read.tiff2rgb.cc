/* 
 Authors
 Martin Schlather, schlather@cu.lu

 library for modeling water flux and solute transport:
   makes TIFFReadRGBAImage available in R
    
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

#ifdef unix

extern "C" {
#include <R.h>
#include <Rdefines.h>
#include <assert.h>
#include <tiffio.h>
#include "SoPhy.h"
}

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
  uint32 Length, Width, total2, total3, i, *lr;
  unsigned long int total;
  int x,Error, *r;
  TIFF* file;
  SEXP raster;
  SEXP dim;

  assert(sizeof(x)==4);
  lr = NULL;

  if ((file = TIFFOpen(STRING_VALUE(filename), "r")) == NULL) {
    Error = 1;
    goto ErrorHandling2;
  } 
  if ((Error=TIFFGetField(file, TIFFTAG_IMAGELENGTH, &Length))!=1) {
    goto ErrorHandling2;
  }
  if ((Error=TIFFGetField(file, TIFFTAG_IMAGEWIDTH, &Width))!=1) {
    goto ErrorHandling2;
  }
  
// printf("width=%d   height=%d \n", Width,Length);
  
  PROTECT(dim=allocVector(INTSXP, 3));
  INTEGER_POINTER(dim)[0] = Width;
  INTEGER_POINTER(dim)[1] = Length;
  INTEGER_POINTER(dim)[2] = 4;
  PROTECT(raster=allocArray(INTSXP, dim));
  r = (int*) INTEGER_POINTER(raster);
  total = Length * Width; 

//printf("total %d %e %e %e %d\n",
// total, (double) total, (double) Width, (double) Length,  sizeof(u_long));

  lr = (uint32*) malloc(total * sizeof(uint32));
  if ((Error = TIFFReadRGBAImage(file, Width, Length, lr, 1))!=1) {
    Error = 999;
    goto ErrorHandling;
  }
  total2 = total * 2;    
  total3 = total * 3;   
  for (i=0; i<total; i++) {
    r[i + total3] = TIFFGetA(lr[i]);
    r[i + total2] = TIFFGetB(lr[i]);
    r[i + total ] = TIFFGetG(lr[i]);
    r[i] = TIFFGetR(lr[i]);
  }
  UNPROTECT(2);

  free(lr);
  return raster;

 ErrorHandling:
  UNPROTECT(2);
 ErrorHandling2:
  PRINTF("error nr %d occured. \n", Error);
  raster=allocVector(INTSXP, 0);
  if (lr != NULL) free(lr);
  return raster;
}

#endif // unix
