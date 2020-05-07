#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/*
extern void PeriodCPT(
  double  *data,       //Periodic time series data
  int     *time,       //Within period time index
  int     *n,          //Length of data
  int     *N,          //Period length
  int     *l,          //Minimum segment length
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


extern void PeriodCPT(
    double  *data,       //Periodic time series data
    int     *time,       //Within period time index
    int     *n,          //Length of data
    int     *N,          //Period length
    int     *l,          //Minimum segment length
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
){

  *err = 0;
  //Initialise draw output to blank character
  for(int i = 0; i < (*nchains * *niter * *Mmax); i++){
    if(i % *Mmax == 0){
      draw[i] = 1;
    }else{
      draw[i] = *blank;
    }
  }

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

