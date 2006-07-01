
/* 
 Authors
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 Collection of system specific auxiliary functions

 Copyright (C) 2001 -- 2006 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

/*
#ifdef WIN32
// #define WIN32_LEAN_AND_MEAN
#define VC_EXTRALEAN
#include <windows.h>
#endif
*/

// achtung! windows.h zusammen mit <Rmath.h oder R.graphics>
// gibt warnung, da ERROR mehrfach definiert !
// deshalb auch in auxiliary.h nicht basic.h einbinden
#include "win_linux_aux.h"
 #include <stdlib.h>
#include <tiffio.h>


int getTIFFinfo(char *filename, 
		long unsigned int *length, long unsigned int *width) {
  TIFF* file;
  uint32 Length, Width;
  int Error;
  if ((file = TIFFOpen(filename, "r")) == NULL) {
    Error = 1;
    goto ErrorHandling1;
  } 
  if ((Error=TIFFGetField(file, TIFFTAG_IMAGELENGTH, &Length))!=1) {
    goto ErrorHandling;
  }
  if ((Error=TIFFGetField(file, TIFFTAG_IMAGEWIDTH, &Width))!=1) {
    goto ErrorHandling;
  }
  TIFFClose(file);
  *length = Length;
  *width = Width;
  return 0;

 ErrorHandling:
  TIFFClose(file);
 ErrorHandling1:
  return Error;
}

int getTIFFimage(char *filename, long unsigned int length,
		 long unsigned int width, int *r) {
  TIFF* file;
  uint32 total2, total3, i, *lr, Width, Length;
  int Error;
  unsigned long int total;


  Length = length;
  Width = width;
  lr = NULL;
  if ((file = TIFFOpen(filename, "r")) == NULL) {
    Error = 1;
    goto ErrorHandling1;
  } 
  total = Length * Width; 
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
  free(lr);
  TIFFClose(file);
  return 0;

 ErrorHandling:
  TIFFClose(file);
 ErrorHandling1:
  if (lr != NULL) free(lr);
  return Error;
}
