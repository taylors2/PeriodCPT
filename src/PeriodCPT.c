#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <string.h>

/*
extern void PeriodCPT(
  double  *data,       //Periodic time series data
  int     *time,       //Within period time index
  int     *n,          //Length of data
  int     *N,          //Period length
  int     *minseglen,  //Minimum segment length
  int     *Mmax,       //maximum number of possible pcpts: floor(N/l)
  char   **Mdist,      //Name of pcpt total prior distribution
  double  *Mhyp,       //Hyper-parameter(s) for pcpt total prior
  int     *spread,     //Hyper-parameter for pcpt location prior
  char   **Pdist,      //Name of sampling distribution
  double  *Phyp,       //Hyper-parameter(s) for segment parameter priors
  int     *inits,      //Chain inital pcpt values
  int     *nchains,    //Number of chains
  int     *nburn,      //Length of burn-in period
  int     *niter,      //Length of chain (after burn-in)
  int     *cache,      //Cache size for the proposal distribution
  int     *quiet,      //Progress bar logical flag
  int     *blank,      //Holding number of blank or NA spaces in inits and draw
  int     *err,        //Error flag
  int     *draw        //Vectorised MCMC output (append each chain)
);
*/

#include "PeriodCPT_LOOKUP.h"

extern void PeriodCPT(
    double  *data,       //Periodic time series data
    int     *time,       //Within period time index
    int     *n,          //Length of data
    int     *N,          //Period length
    int     *minseglen,  //Minimum segment length
    int     *maxM,       //maximum number of possible pcpts: floor(N/l)
    char   **Mdist,      //Name of pcpt total prior distribution
    double  *Mhyp,       //Hyper-parameter(s) for pcpt total prior
    double  *spread,     //Hyper-parameter for pcpt location prior
    char   **Pdist,      //Name of sampling distribution
    double  *Phyp,       //Hyper-parameter(s) for segment parameter priors
    int     *inits,      //Chain inital pcpt values
    int     *nchains,    //Number of chains
    int     *nburn,      //Length of burn-in period
    int     *niter,      //Length of chain (after burn-in)
    int     *cache,      //Cache size for the proposal distribution
    int     *quiet,      //Progress bar logical flag
    int     *blank,      //Holding number of blank or NA spaces in inits and draw
    int     *err,        //Error flag
    int     *draw        //Vectorised MCMC output (append each chain)
){

  *err = 0;

  //MAKE LOOKUP TABLES
  double   *g1;
  segval_t *g2;
  g1 = (double *)   calloc( *maxM ,   sizeof(double)   );
  g2 = (segval_t *) calloc( *N * *N , sizeof(segval_t) );
  MAKE_LOOK_TABLES(data, time, n, N, minseglen, maxM, Mdist, Mhyp,
                   spread, Pdist, Phyp, err, g1, g2);
  if(*err != 0) goto ESC1;


  for(int ichain = 0; ichain < *nchains; ichain++){


  }



  //Initialise draw output to blank character
  for(int i = 0; i < (*nchains * *niter * *maxM); i++){
    if(i % *maxM == 0){
      draw[i] = 1;
    }else{
      draw[i] = *blank;
    }
  }


  ESC1:;
  return;
}


void R_init_PeriodCPT_Cfunctions(DllInfo *info){

  R_CMethodDef cMethods[] = {
    {"PeriodCPT", (DL_FUNC) &PeriodCPT, 20},
    {NULL,NULL,0}
  };

  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

//------------------------
//Error codes:
//------------------------
// 1: Mprior distribution not recognised
// 2: Sampling distribution not recognised
