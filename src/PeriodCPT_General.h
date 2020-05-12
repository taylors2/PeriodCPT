#ifndef FILE_GENERAL
#define FILE_GENERAL

#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <math.h>

#define TRUE 1
#define FALSE 0

typedef struct segval segval_t;
typedef struct MCMCitem MCMCitem_t;
typedef struct MCMCchain chain_t;
typedef struct CacheItem cache_t;

typedef double (Mprior_Ptr       )(int,     int*, double*);
typedef double (Samp_Dist_Ptr    )(double*, int,  double*);
typedef double **(Summary_Stats_Ptr)(double*, int*, int*, int*);


struct segval {
  int     cpt;
  int     seglen;
  double  value;
};

struct MCMCitem{
  int         *tau;
  int          j;
  int          m;
  double       value;
  double       prob;
  MCMCitem_t  *prev;
  MCMCitem_t  *next;
};

struct MCMCchain{
  MCMCitem_t  *first;
  MCMCitem_t  *last;
  int          length;
};

struct CacheItem{
  MCMCitem_t *generator;
  int         count;
  chain_t    *chain;
};





//----------------------------

void Progress(int i, int n, int tk, char **str){
  if(i == 0){
    printf("%s |", *str);
  }

  int diff = floor(tk * (i+1) / n) - floor(tk * (i) / n);
  if(diff != 0){
    for(int itk = 0; itk < diff; itk++){
      printf("=");
    }
  }

  if((i+1) == n){
    printf("|\n");
  }

  return;
}




#endif //FILE_GENERAL
