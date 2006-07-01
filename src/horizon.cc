/* 
 Authors
 Martin Schlather, schlather@cu.lu

 library for modeling water flux and solute transport:
    internal update of the definition of the horizons, 
    especially boundaries between horizons
    
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


#include <R.h>
#include <Rdefines.h>
 
#include "SoPhy.h"
  //#include "auxiliary.h"


/*
  to do:
  reduction of size for horizons. 
  to this end 
    crossing of horizons must be checked!
    Best to do it at the very end
 */

//----------------------------------------------------------------
//functions from RandomFields/auxiliary.cc

void orderdouble(double *d, int *pos, int start, int end END_WITH_GSLRNG)
     /* quicksort algorithm, slightly modified:
        does not sort the data, but d[pos] will be ordered 
	NOTE: pos must have the values 0,1,2,...,start-end !
	(orderdouble is a kind of sorting pos according to
	the variable d)
     */
{
  int randpos, pivot, left, right, pivotpos, swap;
  double Dpivot;

  if( start < end )
  {   
    randpos = start + (int) (UNIFORM_RANDOM * (end-start+1));
    pivot = pos[randpos];
    pos[randpos] = pos[start];
    pos[start] = pivot;
    Dpivot = d[pivot];
   
    pivotpos=start; 
    left = start;
    right=end+1;   

    while (left < right) {      
      while (++left < right && d[pos[left]] < Dpivot) pivotpos++;
      while (--right > left && d[pos[right]] > Dpivot) ;      
      if (left < right) {
        swap=pos[left]; pos[left]=pos[right]; pos[right]=swap;
        pivotpos++;
      }
    }
    pos[start] = pos[pivotpos];
    pos[pivotpos] = pivot;
    orderdouble( d, pos, start, pivotpos-1 END_WITH_RANDOM);
    orderdouble( d, pos, pivotpos + 1, end END_WITH_RANDOM);
  }
}

void quicksortdouble(double *d, int start, int end) {
  int i=start, j=end;
  double h, x=d[(start+end)/2];
  do {    
    while (d[i]<x) i++; 
    while (d[j]>x) j--;
    if (i<=j)
      {
	h=d[i]; d[i]=d[j]; d[j]=h;
	i++; j--;
      }
  } while (i<=j);
  if (start<j) quicksortdouble(d, start, j);
  if (i<end) quicksortdouble(d, i, end);
}

//----------------------------------------------------------------


SEXP XGetHorizons(SEXP h) {
  SEXP x, names;
  
  PROTECT(h);
  PROTECT(names = getAttrib(h, R_NamesSymbol));
  SET_VECTOR_ELT(names, 5, mkChar("abc"));

  PROTECT(x = NEW_NUMERIC(1));
  REAL(x)[0] = 4.5; // or use NUMERIC_POINTER, p.35 bottom
  SET_VECTOR_ELT(h,5,x);
  setAttrib(h, R_NamesSymbol, names);
  UNPROTECT(4);
  return h; 
}


SEXP getListElement(SEXP list, char *str){
  SEXP elmt = R_NilValue, names=getAttrib(list, R_NamesSymbol);
  int i;
  for (i=0; i<length(list); i++)
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

int getListNr(SEXP list, char *str){
  SEXP names=getAttrib(list, R_NamesSymbol);
  int i;
  for (i=0; i<length(list); i++)
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0)  return(i);
  return -1;
}

#define round(x) (floor((x) + 0.5))
void ins_border_point(double x, double y, double *(grid[2]), double step, 
		      int *border, int *actborder, int nborder, int *dim) {
  int intx, inty;

  //  if(*actborder>=nborder) printf("act=%d nb=%d\n", *actborder, nborder);

  assert(*actborder < nborder);
  intx = (int) round((x - grid[0][0])/step);
  inty = (int) round((y - grid[1][0])/step);

  //printf(" ins %d: %f %f %d y=%f %f %d [%f %f %f]\n",
  //	   *actborder, 
  //   x, round((x - grid[0][0])/step), intx, 
  //   y, round((y - grid[1][0])/step), inty,
  //   grid[0][0], grid[1][0], step);

  if ((intx>=0) && (intx<dim[0]) && (inty>=0) && (inty<dim[1])) {
    if ((*actborder==0) ||
	(border[*actborder - 1] != intx) || 
	(border[*actborder + nborder - 1] != inty)) {
      //      if (*actborder>0)
      //	printf("        %d %d != %d  %d != %d\n", *actborder,
      //	       border[*actborder - 1], intx, 
      //	       border[*actborder + nborder - 1], inty); 

      // not all multiple border points voided, but most of them
      border[*actborder] = intx;
      border[*actborder+nborder] = inty;
      //      printf("%d %d (%d, %d)", *actborder, nborder, border[*actborder], 
      //      	     border[*actborder+nborder]);
      //      printf(" ok");
      (*actborder)++;
    }
  }
  //  printf("\n");
}


typedef struct chainelmt {
  double x1, y1, x2, y2;
  chainelmt *next;
} chainelmt;

typedef struct chaintype {
  chainelmt *segment, *last;
  double *interval;
  int n, nint, actint;
} chaintype;

void showchain(char *txt, chaintype *chain) {
  chainelmt *i;
  int j;
  
  return;

  PRINTF("\n %s -Inf=%f  n=%d  nint=%d  act=%d ",txt,
	 chain->n, chain->nint, chain->actint);
  if (chain->interval==NULL) PRINTF("interval=NULL"); else {
    PRINTF("[");
    for (j=0; j< chain->nint*2; j++) PRINTF("%.2f ",chain->interval[j]);
    PRINTF("]");
  }
  i=chain->segment;
  if (i!=NULL) {
    i = i->next;
    while (i!=NULL) { 
      PRINTF(" {%.2f,%.2f; %.2f,%.2f}",i->x1,i->y1,i->x2,i->y2);
      i = i->next;
    }
    PRINTF(" NULL\n"); 
  } else PRINTF("segment is NULL\n");
}

void initchain(chaintype *chain) {
  assert((chain->segment==NULL) && (chain->last==NULL) && 
	 (chain->interval==NULL));
  chain->segment = chain->last = (chainelmt*) malloc(sizeof(chainelmt));
  chain->segment->x2 = chain->segment->x1 = RF_INF;
  chain->segment->next = NULL;
  chain->n = 0;
}

void insertchain(chaintype *chain, double x1, double y1, double x2, double y2) {
  double swap;
  chainelmt *dummy, *i;
  char insert[] = "insert";

  assert(chain->segment!=NULL);

  //showchain("insert",chain);

  if (x1>x2) {
    swap=x1; x1=x2; x2=swap; 
    swap=y1; y1=y2; y2=swap;
  }
  (chain->n)++;

  i = chain->segment;
  dummy = (chainelmt*) malloc(sizeof(chainelmt));
  while ((i->next != NULL) && (i->next->x2 < x2)) i=i->next;
  if (i->next==NULL) {
    i->next = chain->last = dummy;
    i->next->next = NULL;
  } else {
    dummy->next = i->next;
    i->next = dummy;
  }
  i = i->next;
  i->x1 = x1;
  i->y1 = y1;
  i->x2 = x2;
  i->y2 = y2;
  //printf("insert (%f,%f) - (%f,%f)\n", x1, y1, x2, y2);

  showchain(insert, chain);


}

void firstinterval(chaintype *chain, double x, double *u, double *l){
  chainelmt *i, *dummy;
  int j, nint2;

  if (chain->n==0) {
    // can happen if x value is outside projection of polygon to the x-axis
    *u = *l = RF_INF;
    return;
  }

  //  printf("%f %d %d %d \n",x,chain->segment,chain->n,chain->interval);
  assert((chain->segment!=NULL) && (chain->interval==NULL));

 /* update chain, i.e. segments whose right end points are left from x 
  are eliminated */
  // showchain("first",chain);
  i = chain->segment;
  while ((i->next!=NULL) && 
	 (i->next->x2 < x)) // the segments are inserted in the list
      //                       by ordering of endpoints of the segments
      //                       hence, it is ensured that all segments
      //                       whose right end point are left from x 
      //                       are erased 
  { 
    // printf(" deleting [%f] (%f %f) -- (%f,%f) \n",x,i->next->x1,i->next->y1,i->next->x2,i->next->y2);
    dummy=i->next;
    i->next = dummy->next;
    free(dummy);
    chain->n--;
  }

  if (chain->n==0) {
    // can happen if x value is outside projection of polygon
    *u = *l = RF_INF;
    return;
  }
  
  assert(chain->n!=0);

  chain->nint = (chain->n + 1) / 2;
  chain->actint = 0;
  nint2 = chain->nint * 2;
  //printf("======= %d %d %d \n", chain->n, chain->nint, nint2);  
  chain->interval = (double*) malloc(sizeof(double) * nint2);

  //showchain("first",chain);

  i = chain->segment->next;
  for (j=0; j<chain->n; j++) {
    //printf("j= %d \n",j);
    assert(i!=NULL);
    assert(i->x1 <= x);
    assert(x <= i->x2);
    if (i->x1==i->x2) { // and  then == i->x1==x 
      chain->interval[j] = (i->y1 < i->y2) ? i->y1 : i->y2;
      // do not include border line
    } else {
      chain->interval[j] = 
	((i->y2 - i->y1) * x - i->y2 * i->x1 + i->y1 * i->x2) / (i->x2 - i->x1);
    }
    i = i->next;
  }
  
  if (nint2 > chain->n) chain->interval[nint2-1] = RF_INF;
  //for (j=0; j<nint2; j++) printf(" %f ",chain->interval[j]); printf("\n");
  quicksortdouble(chain->interval, 0, nint2-1);

  // printf("check sortdouble\n");
  //for (j=0; j<nint2; j++) printf(" %f ",chain->interval[j]); printf("\n");
  for (j=1; j<nint2; j++) assert(chain->interval[j-1] <=chain->interval[j]); 

  *l = chain->interval[2 * chain->actint];
  *u = chain->interval[2 * chain->actint + 1];
  //printf("exit %f %d %d %d %f %f\n",x,chain->segment,chain->n,chain->interval,*l,*u);
  // showchain("first",chain);
}

void nextinterval(chaintype *chain, double *u, double *l){
  assert(chain->segment != NULL);
  if (++(chain->actint) >= chain->nint) {
    assert(chain->interval != NULL);
    free(chain->interval);
    chain->interval = NULL;
    *u = *l = RF_INF;
    return;
  }
  *l = chain->interval[2 * chain->actint];
  *u = chain->interval[2 * chain->actint + 1];
  //printf("next exit (new act) %d %d %f %f\n",chain->actint,chain->nint,*l,*u);
}

void deletechain(chaintype *chain){
  chainelmt *dummy;
  if (chain->interval!=NULL) {
    free(chain->interval);
    chain->interval = NULL;
  }
  while (chain->segment!=NULL) {
    dummy = chain->segment;
    chain->segment = chain->segment->next;
    free(dummy);
  }
  chain->segment = chain->last = NULL;
}


#define GLOBAL_N 5
#define H_GLOBAL 7
#define I_GRID 0
#define I_IDXRF 2
#define I_N 3
#define I_STEP 4 
#define MAX_HORIZON 10 /* == # of register in RandomFields */
#define H_OPT 7 
#define H_MUST 7
#define OPT_N 4
#define I_CUT 0
#define I_BORDER 2
#define I_IDX 3
#define MUST_N 2
#define I_POINTS 0
#define I_TYPE 1

char h_opt[OPT_N][H_OPT] = {"cut.x",  // area containing horizon or polygon
                            "cut.y",  // range of x and y coordinates
                            "border", // border points
                            "idx"     // index which pixels of the clipping
			    //           area belong to horizon or polygon
                           };
char h_must[MUST_N][H_MUST] = {"points", "type"};
char global_must[GLOBAL_N][H_GLOBAL] = 
{"grid.x", "grid.y", "idx.rf", "n", "step"};
char hH[] = "H";
char hP[] = "P";
char hx[] = "x";
char hy[] = "y";


SEXP GetHorizons(SEXP h, SEXP nr) {
// h: the list of all horizons
// nr[1:2]: range of number of horizons for which the values should be updated 
// nr[2]: is always h$n in standard applications
// nr[1]: is 1, if called by 'calculate.horizons' -- to be on the same side
//              or within xswms2d, except
//        h$n, if a new horizon or polygon is drawn    
  SEXP  elmt, Pts, hi, hinew, namesnew, must_sexp[MUST_N];  
  double *(grid[2]), *(pts[2]), step, *xs, gridMstep, maxgridy;
  char text[100]; 
  int i, j, k, n, dim[2], firstnr, lastnr, npts, given_n, *idx_rf, *idx,
    missing_n,
    protct, *(cut[2]), ncut_tot, intcoord, ncut[2], nborder, actborder, *border, 
    *pos, rf_size, opt_pos[OPT_N];  
  chaintype chain;
  
  chain.segment = chain.last = NULL;
  chain.interval = NULL;
  dim[0] = dim[1] = -1;
  step = -1.0; 
  npts = -1;
  idx_rf = NULL;
       
  PROTECT(h); protct=1;
  for (i=0; i<GLOBAL_N; i++) {
    if ((elmt=getListElement(h, global_must[i]))==R_NilValue) {   
      UNPROTECT(protct);
      sprintf(text, "`%s' not in the list", global_must[i]);
      error(text);
    }
    switch(i){
    case I_GRID :
      case I_GRID+1 :
	assert(isVector(elmt) && isReal(elmt) && (LENGTH(elmt)>1));
	grid[i-I_GRID] = REAL(elmt); 
	dim[i-I_GRID] = LENGTH(elmt);
	break;
    case I_IDXRF :
      if (isVector(elmt) && (LENGTH(elmt)==1) && isNumeric(elmt)
	  &&(ISNAN(REAL(elmt)[0]))) { 
	assert((dim[0]>0) && (dim[1]>0)); // i.e make sure that dim is 
        // initialized -- might be relevant if programme is changed 
	SET_VECTOR_ELT(h, getListNr(h, global_must[i]), 
		       allocMatrix(INTSXP, dim[0], dim[1]));
	elmt=getListElement(h, global_must[i]);
      }
      assert(isMatrix(elmt) && isInteger(elmt));
      idx_rf = INTEGER_POINTER(elmt); // !! Vector
      break;
    case I_N : 
      assert(isVector(elmt) && isInteger(elmt) && (LENGTH(elmt)>0));
      n = INTEGER(elmt)[0]; 
      break;
    case I_STEP : 
      assert(isVector(elmt) && isReal(elmt) && (LENGTH(elmt)>0));
      step = REAL(elmt)[0]; 
      break;
    default : assert(false);
    }
  } 
  assert(step>0 && idx_rf!=NULL);
  rf_size = dim[0] * dim[1]; 
  gridMstep = grid[1][0] - step;
  maxgridy = grid[1][dim[1]-1];

  // horizon 0 always the last, and always(!) calculated, see next but one loop
  if ((firstnr = INTEGER(nr)[0] - 1)==0) {
    firstnr=1;
    for (i=0; i<rf_size; i++) idx_rf[i] = 0;
  }
  lastnr = INTEGER(nr)[1];

  i=firstnr;
  if (i<=lastnr) while (true) { // "for (i=firstnr; i<=lastnr; i++)"
    if (i==lastnr) {
      i=0; // horizon 0 always the last one
    }
 
    // check for list element of h[[i]] that must be given
    hi = VECTOR_ELT(h,i);
    for (j=0; j<MUST_N; j++) {
      if ((elmt = must_sexp[j] = getListElement(hi, h_must[j]))==R_NilValue) { 
	if ((i==0) && (j==I_POINTS)) continue;// hor 0: of course no points that 
        //                                       define any boundary!
	UNPROTECT(protct);
	sprintf(text, "`%s' is not an element of horizon[[%d]]", h_must[j], i);
	error(text);
      }
      switch(j) {
      case I_TYPE : break;
      case I_POINTS : 
	if ((Pts = getListElement(elmt,hx))==R_NilValue)
	  strcpy(text, "no x component");
	else {
	  pts[0] = REAL(Pts);
	  npts = LENGTH(Pts);
	  if ((Pts = getListElement(elmt,hy))==R_NilValue)
	      strcpy(text, "no y component");
	  else {
	    pts[1] = REAL(getListElement(elmt,hy));
	    if (LENGTH(Pts)!=npts) 
	      strcpy(text, "lengths of `x' and `y' in `points' do not match");
	    else break;
	  }
	}
	UNPROTECT(protct);
	sprintf(text, "`%s' in horizon[[%d]]", text, i);
	error(text);
      default : assert(false);
      }
    }
    // assert(npts>0); see below!
    
    // check if "cut.x", "cut.y", "border", "idx" are elements of h[[i]]
    // they are not necessarily given (OPT_N==4)
    given_n = 0;
    for (j=0; j<OPT_N; j++) {
      if ((k=opt_pos[j]=getListNr(hi, h_opt[j]))>=0) {
	given_n++;
	switch (j){
	case I_CUT:
	case I_CUT+1 :
	    // if the optional parameters do not have the right format
	    // the programme ends in a rather unfriendly way
	    //
	    // this only happens if the user changes these optional values
	    // or if there is a programming error 
	  assert(isVector(VECTOR_ELT(hi, k)) && isInteger(VECTOR_ELT(hi, k)));
	  break;
	case I_BORDER :    
	  assert((VECTOR_ELT(hi, k)==R_NilValue) ||
		 (isMatrix(VECTOR_ELT(hi, k)) && isInteger(VECTOR_ELT(hi, k))));
	  break;
	case I_IDX :
	  assert((VECTOR_ELT(hi, k)==R_NilValue) ||
		 (isMatrix(VECTOR_ELT(hi, k)) && isLogical(VECTOR_ELT(hi, k))));
	  break;
	}
      }
    }
    if ((missing_n = OPT_N - given_n) > 0) {
      given_n = LENGTH(hi); // not this resetting and that !! 
      PROTECT(namesnew = allocVector(STRSXP, given_n + missing_n));
      PROTECT(hinew = allocVector(VECSXP, given_n + missing_n));
      for (j=0; j<given_n; j++) {
	SET_VECTOR_ELT(hinew, j, VECTOR_ELT(hi, j));
	SET_STRING_ELT(namesnew, j, STRING_ELT(getAttrib(hi, R_NamesSymbol), 
		       j));
      }
      for (j=0; j<OPT_N; j++) if (opt_pos[j]<0) { // i.e. that variable has 
        //                                           not been given 
	opt_pos[j] = given_n++; // not that the meaning of given_n has changed
	//                         it is now LENGTH(hi) at start
	SET_STRING_ELT(namesnew, opt_pos[j], mkChar(h_opt[j]));
      }
      setAttrib(hinew, R_NamesSymbol, namesnew);
      SET_VECTOR_ELT(h, i, hinew);
      hi = VECTOR_ELT(h, i);
      UNPROTECT(2);
    }

 
    if (i==0) {
      // cut[][] returns matrix indices for R matrices

      // special treatment for first horizon -- treatment must be the last
      // one, using indirectly (?!) results from the other horizons
      //  
      // 19.6.03: the below definitions do not depend on other horizons
      //          but are on the very save side; hence the 'cut' boundaries
      //          can be improved
      for (j=0; j<2; j++) {
	k = opt_pos[I_CUT+j];
	SET_VECTOR_ELT(hi, k, allocVector(INTSXP, 2));
	cut[j] = INTEGER_POINTER(VECTOR_ELT(hi, k));
	cut[j][0] = 1;
	cut[j][1] = dim[j];
      }
      break; // end of loop "for (i=firstnr; i<=lastnr; i++)
    }

    //printf(" cut\n");
   
    for (j=0; j<2; j++) {
      k = opt_pos[I_CUT+j];
      SET_VECTOR_ELT(hi, k, allocVector(INTSXP, 2)); 
      // the previous vector will be erased by garbage collection 
      cut[j] = INTEGER_POINTER(VECTOR_ELT(hi, k));
    }

    assert(npts>0);
    if (strcmp(CHAR(STRING_ELT(must_sexp[I_TYPE], 0)), hH)==0) {
      // it is a horizon definition !
      int intcoord;
      cut[0][0] = 1;      // not optimised yet -- this is worth doing it
      cut[1][0] = dim[1]; // * the lower bound of h[[i]] is optimised below
      //                     * starting value for optimisation
      //                     * conservative value would be 1
      cut[0][1] = dim[0]; // not optimised yet -- worth doing it??
      cut[1][1] = dim[1]; // not optimised yet -- worth doing it??
      for (k=0; k<npts; k++) {
        intcoord = (int) (round((pts[1][k] - grid[1][0])/step)) + 1;//+1;27.12.03
        if (intcoord < cut[1][0]) cut[1][0] = intcoord;
      }
      if (cut[1][0] < 1) cut[1][0]=1;
    } else {
      // polygon
      for (j=0; j<2; j++) {
        cut[j][0] = dim[j];
        cut[j][1] = 1;
        for (k=0; k<npts; k++) {
	  intcoord = 1 + (int) (round((pts[j][k] - grid[j][0]) / step));
	  if (intcoord < cut[j][0]) cut[j][0] = intcoord;
          if (intcoord > cut[j][1]) cut[j][1] = intcoord;
        }
        if (cut[j][0]<1) cut[j][0]=1;
        if (cut[j][1]>dim[j]) cut[j][1]=dim[j];
      }
    }

    // create idx
    ncut[0] = cut[0][1] - cut[0][0] + 1;
    ncut[1] = cut[1][1] - cut[1][0] + 1;
    ncut_tot = ncut[0] * ncut[1];

    k = opt_pos[I_IDX];
///////////////////////
    SET_VECTOR_ELT(hi, k, allocMatrix(LGLSXP, ncut[0], ncut[1]));
    idx = LOGICAL_POINTER(VECTOR_ELT(hi, k));
    for (j=0; j<ncut_tot; j++) idx[j]=false;
      
    // sorting the segments of the horizon/polygon boundary by the
    // x-ccordinate of the left point.  
    xs = (double*) malloc(sizeof(double) * npts);
    pos = (int*) malloc(sizeof(int) * npts);
    int nptsM1;
    nptsM1 = npts-1;
    for (j=0; j<nptsM1; j++) {
      pos[j] = j;  // initialisation will give the ordering
      if (pts[0][j] <= pts[0][j+1]) xs[j] = pts[0][j]; else xs[j] = pts[0][j+1];
    }
    orderdouble(xs, pos, 0, nptsM1-1);
    initchain(&chain);
    // just to be sure that the ordering worded
    // 19.6.03 should be deleted later on...  
    for (j=1; j<nptsM1; j++) assert(xs[pos[j-1]] <= xs[pos[j]]);
    
    double x, u, l;
    int iy, jy, upper;
    k = 0;
    x = (double) (cut[0][0] - 1) * step + grid[0][0];
    for (j=cut[0][0]; j<=cut[0][1]; j++, x += step){

      // 27.12.03 : keep following comments:
      //printf("upd %d %d %f %d %d\n",j,cut[0][1],x, k);
      // if (k<nptsM1) printf("** k= %d  %f %f %d \n",k, x,xs[pos[k]],pos[k]);

      // if by the next pixels a new segment is crossed,
      // insert this segment (or these segments) into the list
      // the segments that are left behind are automatically be deleted
      // from the list by 'firstinterval'
      // think of all relevant segments be piled up (with some interspace)
      // mathematics say that now the parts of the y-segments between the 
      // piled boundary segment of the clipping area at the given x-coordinate 
      // alternate between belonging and not belonging to the horizon/polygon  
      while (k < nptsM1 && xs[pos[k]] < x) {
	// xs[pos[k]] < x : left open right closed interval, 
	// cf. 'firstinterval' where a segment is deleted from list only if
	// the right end point is truly left from the current position
	// (( if both sides of the interval were closed and the x-coordinate
	// of a end/start point of a segment equals the current position  
	// than the old segment would not have delete and the new one is
        // already inserted; then the algorithm behaves as if there were
	// two independent segments with interspace 0 in y-direction
        // Exception: for the very first point the interval is also left closed,
	//            see 'firstinterval'
 
	// the selection of the pts works in a trivial way: what ever the value
	// of pos[k] is, this pts[][pts[k]] and the following one is taken.
	// since pos[k] runs through 0...nptsM1-1 all segments are considered
	// (since npts is the number of points, ntpsM1 the number of pairs
	// and the index runs starts at 0)
	insertchain(&chain, 
		    pts[0][pos[k]], pts[1][pos[k]],
		    pts[0][pos[k] + 1], pts[1][pos[k] + 1]);
	k++;
      }
      firstinterval(&chain, x, &u, &l);
      if (false)  // debugging only
	PRINTF("\n j=%d, u=%f l=%f \n",j,u,l);
      while (l < RF_INF) {
	// -1 since cut takes index values for R
	// (1+...) since only interior of region should be set
	///printf("l=%f %f %f\n",l, maxgridy );
	if (l <= maxgridy) {
	  //printf("ok/\n",u);
	  if (u>maxgridy) u = maxgridy;
	  if (l<gridMstep) l = gridMstep;
	  iy = j + (1 + (int) ((l - grid[1][0]) / step)) * dim[0] - 1;
	  jy = j - cut[0][0] +
	    (1 + (int) ((l - grid[1][0]) / step) - (cut[1][0]-1)) * ncut[0];
	  upper = j + ((int) ((u - grid[1][0]) / step)) * dim[0] - 1;
	  if (false)  // debugging only
	  {
	    PRINTF("jy=%d j=%d cut00=%d %d cut10=%d ncut0=%d upper=%d ",
		   jy, j, cut[0][0], (int) ((l - grid[1][0]) / step), cut[1][0],
		   ncut[0], upper);
	    PRINTF("l=%5.2f u=%5.2f\n cut=%d %d %d %d\n",
		   l,u, ncut_tot, ncut[1], cut[1][0],cut[1][1]);
	  }
	  for (; iy <= upper; iy += dim[0], jy+=ncut[0]) {
	    assert((jy>=0) && (jy<ncut_tot));
	    assert((iy>=0) && (iy<rf_size));
	    idx[jy] = true;
	    idx_rf[iy] = i;
	  }
	}
	nextinterval(&chain, &u, &l);
      } 
    }

    // rough number of entries into $border
    nborder = 0;
    for (j=1; j<npts; j++) {
      nborder += 3 + (int) (hypot(pts[0][j]-pts[0][j-1],
				   pts[1][j]-pts[1][j-1]) / step);
    }
    border = (int*) malloc(sizeof(int) * 2 * nborder);

    // collect border points
    actborder = 0;
    for (j=1; j<npts; j++) {
      double delta, deltax, deltay, x, y;
      int int_delta;
      if (pts[0][j]==pts[0][j-1]) {
	deltax = 0.0;
	deltay = step;
	if (pts[1][j] < pts[1][j-1]) deltay = -step;
      } else {
	deltax = 1;
	deltay = (pts[1][j] - pts[1][j-1]) / (pts[0][j] - pts[0][j-1]);
	delta = step / hypot(deltax, deltay);
	if (pts[0][j]<pts[0][j-1]) delta = -delta;
	deltax *= delta;
	deltay *= delta;
      }
      int_delta = 
	(int) (hypot(pts[0][j]-pts[0][j-1], pts[1][j]-pts[1][j-1]) / step);
      x = pts[0][j-1];
      y = pts[1][j-1];
      for (k=0;  k<=int_delta;  k++, x += deltax, y += deltay) {
	ins_border_point(x, y, grid, step, border, &actborder, nborder, dim);
      }
    }
    ins_border_point(pts[0][npts-1], pts[1][npts-1], grid, step, 
		     border, &actborder, nborder, dim);
    //printf("\n[%f %f]\n",pts[0][npts-1], pts[1][npts-1]);

    //for(k=0; k<actborder; k++) 
    //  printf("(%d,%d) ",border[k],border[k + nborder]); printf("\n");

    //   printf("transfer \n");
    // transfer border points to `h'
    k = opt_pos[I_BORDER];
    SET_VECTOR_ELT(hi, k, allocMatrix(INTSXP, actborder, 2));
    for (j=0; j<actborder; j++) {
      int ii;
      //printf("transferring %d %d %d %d\n",k,j,border[j], border[j + nborder]);
      INTEGER_POINTER(VECTOR_ELT(hi, k))[j] = border[j] + 1;
      INTEGER_POINTER(VECTOR_ELT(hi, k))[j + actborder] = border[j + nborder] +1;
      ii = border[j] + border[j + nborder] * dim[0];
      assert((ii>=0) && (ii<rf_size));
      // idx_rf[ii] =  -i;  /// 
      idx_rf[ii] =  i; 
    
      //printf(" %d %d %d %d %d %d\n",
      //   border[j]- cut[0][0] + (border[j + nborder] - cut[1][0]) * ncut[0],
      //    border[j], cut[0][0], border[j + nborder], cut[1][0], ncut[0]);

      // if (border[j + nborder] - (cut[1][0]-1)<0) printf("!!!!");

      ii = border[j] - (cut[0][0]-1) +
	(border[j + nborder] - (cut[1][0]-1)) * ncut[0];
      // printf("%d %d %d    x=%d in [%d %d] y=%d in [%d %d] %d %d\n",
      //     ii, ncut_tot, j, border[j], cut[0][0]-1, cut[0][1]-1,
      //     border[j + nborder],cut[1][0]-1, cut[1][1]-1, nborder,  ncut[0]);
      assert((ii>=0) && (ii<ncut_tot)); 
      idx[ii] = false;
    }
    free(xs);
    free(border);
    free(pos);
    deletechain(&chain);
    i++; // for (i=firstnr; i<=lastnr; i++) {
  } // i

  // h idx, cutx, cuty updaten
  // check ueber alle n definitionen
  UNPROTECT(protct);
  return h; 
}

void preproot(double *mat, int *cond, double *uptake, int *lmat,
	      double *beta, int *pl, 
	      int *lpl,  // number of root pixels counted for each plant
	      // note that lpl can be greater than np since the same grid point
	      // can be occupied by roots of different plants. If the
	      // condition types are different, one main condition is determined
	      // since only one type of condition is allowed.
	      // Further, for different roots of the same type at the same grid
	      // point, the value of the condition has to be determined.
	      int *np, // number of possible grid points (length.x * length.y)
	      int *sequ,
	      int *diridx, int *neuidx, int *atmidx, 
	      double *dirval, double *neuval, double *atmval,
	      double *root, int *atmadd, int *Error) {
    // dir=dirichlet, neu=neumann, atm=atmospheric
    // diridx, neuidx and atmidx : indicator which condition should be used
    // dirval, neuval, atmval : value which is used if the condition is indicated
    // root: for each point the values P0, P2H, P2L, P3, r3H, r2L which defines
    //       the trapezoid curve

#define NR_COND 4 // the number of conditions we have to consider: 
    //               dir, neu, atm, none
#define NONE 3 // index for none
    
    int *Z, // keeps track how many conditions are found of which kind
	//  in each node 
      condition, roottype, mat_segment, *(idx[NR_COND]), i, j, k, endfor,
      max, root_segment, coord;
  
  *Error = 0;
  // idx is an array of vectors
  idx[0] = diridx; // dirichlet
  idx[1] = neuidx; // neumann
  idx[2] = atmidx; // atmospherical
  
  if ((Z = (int*) calloc(sizeof(int), *np * NR_COND)) == NULL) {
    *Error = 1; goto ErrorHandling; }
  endfor = *np * 6;
  for (i=0; i<endfor; ) {
    root[i++] = RF_NEGINF; // P0, the higher the better for the plants
    root[i++] = RF_INF;    // P2H, the smaller the better for the plants
    root[i++] = RF_INF;    // P2L, the smaller the better for the plants
    root[i++] = RF_INF;    // P3, the smaller the better for the plants
    i += 2;  // r2H, r2L
  }
  for (i=0; i<*np; i++) {
    neuval[i] = atmval[i] = 0;
    dirval[i] = RF_INF;
  }
  for (i=0; i<*lpl; i++) {
    roottype = pl[i + *lpl] - 1;
    assert((roottype >= 0) && (roottype < *lmat));
    condition = cond[roottype] -1;//condition is 0 for dirichlet, etc. see below 
    coord = pl[i]; // although matrix, elements are numbered as if vector 
    assert((coord>=0) && (coord<*np));
    assert((condition >= 0) && (condition <= 3));
    assert(coord * NR_COND + condition < *np * NR_COND);
    Z[coord * NR_COND + condition]++;
    switch(condition) {
    case 0 : // dirichlet
      if (dirval[coord] > uptake[roottype]) 
	dirval[coord] = uptake[roottype];
      break;
    case 1 : // neumann
      neuval[coord] += uptake[roottype];
      break;
    case 2 : // atmospheric
      mat_segment = roottype * 6;
      root_segment = coord * 6;
      if (root[root_segment] < mat[mat_segment])
	root[root_segment] = mat[mat_segment];
      if (root[root_segment + 3] > mat[mat_segment + 3])
	root[root_segment + 3] = mat[mat_segment + 3];
      for (j=1; j<=2; j++) {
	if (root[root_segment + j] > mat[mat_segment + j]) {
	  root[root_segment + j] = mat[mat_segment + j];
	  root[root_segment + j + 3] = mat[mat_segment + j + 3];
	}
      }
      if (*atmadd) atmval[coord] += beta[i]; else 
	if (atmval[coord]<beta[i]) atmval[coord] = beta[i];
      break;
    case 3 : // none 
      break; 
    default :
      assert(false);
    }
  }


  for (j=i=0; i<*np; i++, j+=NR_COND) {
    max=0;
    assert(j + NR_COND - 1 < *np * NR_COND);
    for (k=0; k<NR_COND; k++) if (Z[j + k] > max && k!=NONE) max = Z[j + k];
    if (max > 0) {
      for (k=0; k<NR_COND; k++) {
	if (Z[j + sequ[k]] == max && k!=NONE) {
	  assert((sequ[k] >= 0) && (sequ[k]<NR_COND));
	  idx[sequ[k]][i] = true;
	  break; // ensures that only one idx[] get the value true.
	}
      }
    }    
  }
  
  free(Z);
  return;
  
 ErrorHandling:
  if (Z != NULL) free(Z);
  return;
}






