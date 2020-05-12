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
#include "PeriodCPT_MCMC.h"
#include "PeriodCPT_MCMCgeneric.h"
#include "PeriodCPT_Cache.h"

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
//  int NTICKS = 20;

  //MAKE LOOKUP TABLES
  double   *g1;
  segval_t *g2;
  g1 = (double *)   calloc( *maxM ,   sizeof(double)   );
  g2 = (segval_t *) calloc( *N * *N , sizeof(segval_t) );
  MAKE_LOOK_TABLES(data, time, n, N, minseglen, maxM, Mdist, Mhyp,
                   spread, Pdist, Phyp, err, g1, g2);
  if(*err != 0) goto ESC1;
/*
  int save;
  char str[50];
  char *progress_text = str;
  for(int ichain = 0; ichain < *nchains; ichain++){
    //initialise chain
    MCMCitem_t *current;
    current = GetInital_tau(&(inits[ichain * *maxM]), N, maxM, blank, g1, g2);
    cache_t **PROPCACHE;
    PROPCACHE = Make_Cache_List(*cache);
    int ncache = 0;

    //Run burn
    save = FALSE;
    sprintf(progress_text, "Chain %d/%d (burn-in)  : ",ichain, *nchains);
    PeriodCPT_MCMC(current, nburn, save, &(draw[ichain * *maxM * *niter]),
      N, maxM, minseglen, g1, g2, PROPCACHE, &ncache, cache,
      quiet, &progress_text, blank, &NTICKS, err);

    //Run chain, syncro exporting of sample to output
    save = TRUE;
    sprintf(progress_text, "Chain %d/%d (iteration): ",ichain, *nchains);
    PeriodCPT_MCMC(current, niter, save, &(draw[ichain * *maxM * *niter]),
      N, maxM, minseglen, g1, g2, PROPCACHE, &ncache, cache,
      quiet, &progress_text,blank, &nticks, err);

    //Clean-up
    Delete_Cache_List(PROPCACHE, *cache);
    Delete_MCMCitem(current);

  }

  //If requested, calculate mode estimates and summaries
  //  Must return void and be able to call from R
*/


  ESC1:;
  free(g2);
  free(g1);

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
