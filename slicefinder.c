#include <stdio.h>
#include <math.h>
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


int input(const char *fn, unsigned **xp, int *np);
int pathlen(unsigned **D, unsigned nd, unsigned i, unsigned j);
int printpth(unsigned **D, unsigned *X, unsigned nd, unsigned ii, unsigned jj);
int removeslice(unsigned **D, unsigned *S, unsigned nd, unsigned ns);

sequence *initseq(unsigned *X, unsigned np);
void printseq(sequence *Q);
slice *initslice(sequence *Q, unsigned ii, unsigned jj, unsigned len);
void printslice(sequence *Q, slice *slc);



int main(int argc, char **argv)  {
  int i, j, np = 0;
  unsigned *X = NULL;

  // read sequence from file
  input("invals.txt", &X, &np);
  sequence *Q = initseq(X, np);
  printseq(Q);
  slice *slp = NULL;

  // all path lengths
  int nd = Q->nd;
  for(i=0;i<nd;i++)  {
    for(j=0;j<nd-i;j++)  {
      unsigned pl = 1 + pathlen(Q->D, nd, i, j);   // number of intervals
      if (pl >= 2)  {
        printf("----> len = %d : ", pl+1);
        printpth(Q->D, Q->X, nd, i, j);
        slp = initslice(Q, i, j, pl);
        printslice(Q, slp);
      }
    }
  }
}



sequence *initseq(unsigned *X, unsigned np)  {
  sequence *Q = (sequence *)malloc(sizeof(sequence));
  Q->X = X;
  Q->np = np;
  Q->nd = np-1;
  int i, j;


  int nd = Q->nd;    // convenience

  // D[i][j] is distance from X[i] to X[i+(j+1)]
  unsigned **D = (unsigned **)malloc(nd*sizeof(unsigned *));
  for(i=0;i<nd;i++)  {
    D[i] = (unsigned *)malloc(nd*sizeof(unsigned));     // really only need i
  }

  // compute diffs
  for(i=0;i<nd;i++)  
    D[i][0] = X[i+1] - X[i];

  for(i=1;i<nd;i++)  {
    for(j=0;j<nd-i;j++)  
      D[j][i] = D[j+1][i-1] + D[j][0];
    for(j=nd-i;j<nd;j++) 
      D[j][i] = 0;
  }

  Q->D = D;
  Q->S = NULL;
  return Q;
}

void printseq(sequence *Q)  {

  unsigned **D = Q->D;
  unsigned nd = Q->nd;
  unsigned np = Q->np;

  printf("Sequence of length %u\n\n", np);

  // view sequence
  for(int i=0;i<np;i++)  {
    printf("%u ", Q->X[i]);
  }
  printf("\n\n");

  // view diffs
  for(int i=0;i<nd;i++)  {
    for(int j=0;j<nd-i;j++)  {
      printf("%10u", D[i][j]);
    }
    printf("\n");
  }

}


// Read sequence from file. Eventually, initialize sequence from a SEXP.
int input(const char *fn, unsigned **xp, int *n)  {

  FILE *fp = fopen(fn, "r");
  if (fp == NULL)  {
    printf("error: can't open file %s\n", fn);
    exit(1);
  }

  unsigned nn = 0;
  int nnn;
  fscanf(fp, "%d", &nnn);

  nn = (unsigned) nnn;
  unsigned *X = (unsigned *)malloc(nn*sizeof(unsigned));

  for(int i=0;i<nn;i++)  {
    fscanf(fp, "%d", &nnn);
    X[i] = (unsigned) nnn;
  }
  fclose(fp);
  *n = nn;
  *xp = X;

  return 0;

}

// The sequence is easy to compute from the indices, but the 
// indices are not easy to compute from the points.
// I.e., if we have S = {i1, i2, i3, ... } then it's 
// easy to get the slice points = {X[S[i1]], X[S[i2]], ... }
// But from the points {X1, X2, ... } it's necessary
// to do a search to find i1, i2, etc. 

// Return length of slice from X[i] with step D[i][j] 
// D[i][j] is the difference between X[i] and X[i+j+1]

int pathlen(unsigned **D, unsigned nd, unsigned i, unsigned j)  {
  unsigned val = D[i][j];

  if ((i + j + 1) > nd) 
    return 0;

  // special case: last point of sequence is in slice
  if ((i + j + 1) == nd)     
    return 1;

  // if the step size (val) is in row (i+j+1)
  // then the length of the slice is one more
  // then the length of the rest of the slice.

  for(int k=0;k<(nd - (i+j+1));k++)  {
    if (D[i+j+1][k] == val)  {
      return (1 + pathlen(D, nd, i+j+1, k));
    }
  }
  return 0;

}

// private function for debugging only
int printpth(unsigned **D, unsigned *X, unsigned nd, unsigned ii, unsigned jj)  {

  unsigned i = ii, j = jj;
  unsigned val = D[i][j];
  unsigned found = 0;
  printf("%4u ", X[i]);
  i = (i + j + 1);

  while (1)  {
    if (i > nd)  {
      printf("warning: this should not happen\n");
      return 0;
    }

    if (i == nd)  { 
      printf("%4u ", X[nd]);
    }

    // find column number for diff val in row i
    found = 0;
    for(int k=0;found == 0 && k<nd-i;k++)  {
      if (D[i][k] == val)  {
        printf("%4u ", X[i]);
        i = i + k + 1;
        found = 1;
      }
    }

    if (!found)  {
      printf("%4u ", X[i]);
      break;
    }

  }
  printf("\n");
  return 0;

}

// Find the slice {start = X[i], step = (X[i+j]-X[i])}
// Need to use pathlen first to find the length to allocate
slice *initslice(sequence *Q, unsigned ii, unsigned jj, unsigned len)  {

  unsigned **D = Q->D;
  int nd = Q->nd;
  unsigned step = D[ii][jj];
  unsigned found = 0;
  unsigned i = ii, j = jj;

  // len is the number of steps; number of points is one more than that.
  slice *slc = (slice *)malloc(sizeof(slice));
  unsigned *S = (unsigned *)malloc((len+1)*sizeof(unsigned));

  S[0] = i;
  int s = 1;
  i = (i + j + 1);

  while (1)  {
    if (i > nd)  { 
      printf("warning: this should not happen\n");
      return 0;
    }

    // special case: last point in slice
    if (i == nd)  {
      S[s++] = nd;
    }

    found = 0;
    for(int k=0;found == 0 && k<nd-i;k++)  {
      if (D[i][k] == step)  {
        S[s++] = i;
        i = i + k + 1;
        found = 1;
      }
    }

    if (!found)  {
      S[s++] = i;
      break;
    }
  }

  slc->start = Q->X[S[0]];
  slc->step = step;
  slc->stop = Q->X[S[s-1]];
  slc->ns = s;
  slc->S = S;

  if (s != len+1)  {
    printf("error: found %u points but len is %u:  ", s, len+1);
    printf("start:stop:step = %u:%u:%u\n", slc->start, slc->stop, slc->step);
  }

  return slc;
}

void printslice(sequence *Q, slice *slc)  {
  printf("start:stop:step = %u:%u:%u\n", slc->start, slc->stop, slc->step);
  for(int i=0;i<slc->ns;i++)  
    printf("%u ", slc->S[i]);
  printf("  [");
  for(int i=0;i<slc->ns;i++)  
    printf("%u ", Q->X[slc->S[i]]);
  printf("]\n\n");
}


// Remove slice S of length ns from D0 of size n0 and return D1 (of size n0-ns)
//   S is an array of *indices*, e.g., S = [2, 4, 6] means remove rows [2, 4, 6] and
//       fix up the other rows to account for the removal 
//   n0 is the number of rows in D0, which is one less than the size of the current
//       sequence.
// Assume: all S are less than nd

int removeslice(unsigned **D, unsigned *S, unsigned nd, unsigned ns)  {

  // remove diagonally
  for(int i=0;i<ns;i++)  {
    int s = S[i];
    for(int l=0;l<s;l++)  { 
      D[l][s-l-1] = 0;
    }
  }

  // squeeze out obsolete rows
  int m = 0;    // pull up this many rows
  int j = 0;    // index into S

  for(int i=0;i<nd;i++)  {
    if (i+m == S[j])  {
      j++;     // look for next j
      m++;     // pull back one more row
    }
    if (m == 0)
      continue;

    // updated row i is original row i+m
    if (i+m >= nd)  {
      for(int k=0;k<nd-i;k++)  {
        D[i][k] = 0;
      }
    }
    else  {
      for(int k=0;k<nd-(i+m);k++)  {
        D[i][k] = D[i+m][k];
      }
      for(int k=nd-(i+m);k<nd-i;k++)  {
        D[i][k] = 0;
      }
    }
  }

  // free last rows
  for(int i=0;i<ns;i++)  {
    printf("free row %u\n", nd-i-1);
    free(D[nd-i-1]);
  }

  // squeeze out 0's
  for(int i=0;i<nd-ns;i++)  {
    int m = 0;
    for(int k=0;k<nd-ns-i;k++)  {
      if (D[i][k] == 0)  {
        while(D[i][k+m] == 0)
          m++;
        D[i][k] = D[i][k+m];
        if (m > 0)
          D[i][k+m] = 0;
      }
    }
  }

  // DEBUG: (remove) check for non-zeroed out elements
  for(int i=0;i<nd-ns;i++)  {
    for(int j=nd-ns;j<nd;j++)  {
      if (D[i][j] != 0)  
        printf("Non-zero (%u,%u)\n", i, j);
    }
  }

  return 0;
}

