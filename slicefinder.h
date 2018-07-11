#ifndef __SLICEFINDER__
#define __SLICEFINDER__

#include <stdio.h>
#include <stdlib.h>

typedef struct _slice  {
  unsigned start;
  unsigned stop;
  unsigned step;
  unsigned ns;
  unsigned *S;    // indices not values
} slice;

typedef struct _sequence  {
  unsigned **D;
  unsigned *X;
  unsigned np;    // npoints
  unsigned nd;   // nintervals
  slice *S;
} sequence;

// auxiliary routine for testing, read sequence from file.
// replace with routine that obtains sequence from a SEXP.
int input(const char *fn, unsigned **xp, int *np);

sequence *initseq(unsigned *X, unsigned np);
void printseq(sequence *Q);
slice *initslice(sequence *Q, unsigned ii, unsigned jj, unsigned len);
void printslice(sequence *Q, slice *slc);

int allslices_start_given(sequence *Q, slice ***R, unsigned *nr, unsigned start);

// !!! not fully implemented yet !!!
int removeslice(sequence *Q, slice *slc);  


#endif
