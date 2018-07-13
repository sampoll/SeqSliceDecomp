#include <stdio.h>
#include <stdlib.h>
#include "slicefinder.h"

void check(unsigned *Y, slice **slices, unsigned ns, unsigned *rpoints, unsigned nr, unsigned np);

int main(int argc, char **argv)  {
  unsigned *X = NULL;
  int np = 0;

  int nrpoints = -1, nslices = -1; 
  unsigned *rpoints = NULL;
  slice **slices = NULL;

  if (argc < 2)  {
    printf("usage: slicefinder <infile>\n");
    exit(1);
  }
  input(argv[1], &X, &np);

  // make a copy of X to check afterwards
  unsigned *Y = (unsigned *)malloc(np*sizeof(unsigned));
  for(int i=0;i<np;i++)  {
    Y[i] = X[i];
  }

  // slice decomposition
  algorithm1(X, np, &rpoints, &nrpoints, &slices, &nslices);

  /*
  // Note: indices in slices[] are nonsensical because they
  // refer to different intermediate sequences. 
  printf("Slices {");
  for(int i=0;i<nslices;i++)  {
    printf(" %u:%u:%u   ", slices[i]->start, slices[i]->stop, slices[i]->step);
  }
  printf("}\nRemaining points to put in trivial slices:  ");
  for(int i=0;i<nrpoints;i++)  {
    printf("%u ", rpoints[i]);
  }
  printf("\n");
  */

  check(Y, slices, nslices, rpoints, nrpoints, np);

  int ntotal = nslices + nrpoints / 2;
  if (nrpoints % 2)
    ntotal++;

  char **sres = (char **)malloc(sizeof(ntotal));
  int j = 0;
  for(int i=0;i<nslices;i++)  {
    unsigned sz = snprintf(NULL, 0, "%u:%u:%u", slices[i]->start, slices[i]->stop, slices[i]->step);
    sres[j] = (char *)malloc((sz+1)*sizeof(char));
    sprintf(sres[j++], "%u:%u:%u", slices[i]->start, slices[i]->stop, slices[i]->step);
  }
  for(int i=0;i<nrpoints/2;i++)  {
    unsigned start = rpoints[2*i];
    unsigned stop = rpoints[2*i+ 1];
    unsigned step = stop-start;

    unsigned sz = snprintf(NULL, 0, "%u:%u:%u", start, stop, step);
    sres[j] = (char *)malloc((sz+1)*sizeof(char));
    sprintf(sres[j++], "%u:%u:%u", start, stop, step);
  }
  if (nrpoints % 2)  {
    unsigned start = rpoints[nrpoints-1];
    unsigned stop = start;
    unsigned step = 1;

    unsigned sz = snprintf(NULL, 0, "%u:%u:%u", start, stop, step);
    sres[j] = (char *)malloc((sz+1)*sizeof(char));
    sprintf(sres[j++], "%u:%u:%u", start, stop, step);
  }
  printf("Sequence size = %u\nTotal slices = %u\n", np, j);
  for(int i=0;i<ntotal;i++)  {
    printf("%s ", sres[i]);
  }

  for(int i=0;i<nslices;i++)  {
    free(slices[i]->S);
    free(slices[i]);
  }
  free(slices);
  free(rpoints);
  free(Y);

}


int compare(const void *a, const void *b)  {  return ( *(unsigned *)a - *(unsigned *)b );  }
// Check that no points are lost and no spurious points gained
void check(unsigned *Y, slice **slices, unsigned ns, unsigned *rpoints, unsigned nr, unsigned np)  {
  unsigned *Z = (unsigned *)malloc(np*sizeof(unsigned));
  int nz = 0;
  unsigned toomany = 0;

  for(int i=0;i<ns;i++)  {
    for(int s=slices[i]->start;s<=slices[i]->stop;s += slices[i]->step)  {
      if (nz < np)  
        Z[nz++] = s;
      else  
        toomany = 1;
    }
  }
  for(int i=0;i<nr;i++)  {
      if (nz < np)  
        Z[nz++] = rpoints[i];
      else  
        toomany = 1;
  }

  qsort(Z, nz, sizeof(unsigned), compare);

  int ok = 1;

  if (toomany)  {
    printf("error: too many points\n");
    ok = 0;
  }
  else if (nz != np) {
    printf("error: too few points\n");
    ok = 0;
  }
  else  {
    for(int i=0;i<np;i++)  {
      if (Y[i] != Z[i])
        ok = 0;
    }
  }

  if (!ok)  {
    for(int i=0;i<np;i++)  {
      printf("%u ", Y[i]);
    }
    printf("\n");
    for(int i=0;i<nz;i++)  {
      printf("%u ", Z[i]);
    }
    printf("\n");
  }
  else  {
    printf("Solution OK\n");
  }

  free(Z);

}

// Note: when a slice is removed, the 
// complete set of slices has to be recomputed.
// Because, for example, if we remove 2:8:2, 
// then 3:9:3 is invalidated.

// Is it faster to traverse the set of 
// slices and split up all slices that
// have been invalidated? Or is it 
// faster to simply recompute all slices?

// The former seems complicated because, e.g.,
// removing {18, 24, 30, 36, 42, 48, 54} 
// invalidates {16, 24, 32, 40, 48, 56} and
// leaves {16, 32, 40, 56}, which is not a slice

// For now, recompute all slices after every 
// slice removal. The heuristic is that it is
// hoped that removing longest slices 
// reduces the set so quickly that the 
// recomputation becomes smaller fast.

