#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <string.h>

#include "PeriodCPT_LOOKUP.h"
#include "PeriodCPT_MCMC.h"
#include "PeriodCPT_MCMCgeneric.h"
#include "PeriodCPT_Cache.h"

//---------------------------------------------------------

extern void PeriodCPT_RJMCMC(
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

  //Rprintf("PeriodCPT\n");
  *err = 0;
  int NTICKS = 20;

  //MAKE LOOKUP TABLES
  double   *g1;
  segval_t *g2;
  g1 = (double *)   my_calloc( *maxM ,   sizeof(double)   );
  g2 = (segval_t *) my_calloc( *N * *N , sizeof(segval_t) );
  MAKE_LOOK_TABLES(data, time, n, N, minseglen, maxM, Mdist, Mhyp,
                   spread, Pdist, Phyp, err, g1, g2);
  //PrintLookup(g1, g2, maxM, N);

  if(*err != 0) goto ESC1;

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
    sprintf(progress_text, "Chain %d/%d (burn-in)  : ",ichain+1, *nchains);
    PeriodCPT_MCMC(current, nburn, save, &(draw[ichain * *maxM * *niter]),
                   N, maxM, minseglen, g1, g2, PROPCACHE, &ncache, cache,
                   quiet, &progress_text, blank, &NTICKS, err);

    //Run chain, syncro exporting of sample to output
    save = TRUE;
    sprintf(progress_text, "Chain %d/%d (iteration): ",ichain+1, *nchains);
    PeriodCPT_MCMC(current, niter, save, &(draw[ichain * *maxM * *niter]),
                   N, maxM, minseglen, g1, g2, PROPCACHE, &ncache, cache,
                   quiet, &progress_text,blank, &NTICKS, err);

    //Clean-up
    Delete_Cache_List(PROPCACHE, *cache);
    Delete_MCMCitem(current);

  }

  //If requested, calculate mode estimates and summaries
  //  Must return void and be able to call from R



  ESC1:;
  my_free(g2);
  my_free(g1);

  return;
}

//------------------------------------------------------------


extern void PeriodCPT_EvaluateSufficientStats(
    double  *data,       //Periodic time series data
    int     *time,       //Within period time index
    int     *n,          //Length of data
    int     *N,          //Period length
    char   **Mdist,      //Name of pcpt total prior distribution
    char   **Pdist,      //Name of sampling distribution
    double  *Phyp,       //Hyper-parameter(s) for segment parameter priors
    int     *draws,      //m freq tau1, ...
    int     *ncdraws,    //length of draw per case
    int     *nrdraws,    //number of cases in draw
    int     *nSuffStats, //number of (posterior) sufficient statistics
    double  *blank,      //Holding number of blank or NA spaces in inits and draw
    int     *err,        //Error flag
    double  *out         //Vectorised output
){

  *err = 0;

  //Grab distribution related functions
  Mprior_Ptr *Mprior = NULL;
  Samp_Dist_Ptr *Samp_Dist = NULL;
  Summary_Stats_Ptr *Summary_Stats = NULL;
  Sufficient_Stats_Ptr *Sufficient_Stats = NULL;
  Get_Functions(Mdist, Pdist, &Mprior, &Samp_Dist, &Summary_Stats, err);
  Get_SuffStats_Function(Pdist, &Sufficient_Stats, err);
  if(*err != 0){
    goto ESC1;
  }

  //flood output with blank character
  for(int i = 0; i < *nSuffStats * (*ncdraws-2) * (*nrdraws); i++){
    out[i] = *blank;
  }


  double **SumStats = Summary_Stats(data, time, n, N);
  int numSumStats = (int)SumStats[0][0];
  /*
  Rprintf("nSumStats: %d\n",numSumStats);
  for(int i = 0; i < *N; i++){
    Rprintf("  %d) SumStats - ",i);
    for(int k = 0; k<numSumStats;k++){
      Rprintf("%d: %f ",k+1, SumStats[k+1][i]);
    }
    Rprintf("\n");
  }

  int m, *tau, thistau, prevtau;
  for(int i = 0; i < *nrdraws; i++){
    m = draws[i * *ncdraws];
    tau = (int*)my_calloc(m, sizeof(int));
    for(int k = 0; k < m; k++){
      tau[k] = draws[i * *ncdraws + 2 + k];
    }

    for(int seg = (m-1); seg >= 0; seg--){
      thistau = tau[seg];
      if(seg == 0){
        prevtau = tau[m-1];
      }else{
        prevtau = tau[seg - 1];
      }
    }
    Rprintf("%d] m: %d tau: (%d",i,m,tau[0]);
    for(int k=1; k<m; k++) Rprintf(", %d",tau[k]);
    Rprintf(")\n");

    my_free(tau);
  }*/


  int m, *tau, thistau, prevtau;
  //int freq;
  double *thisStats, *thisSuff;
  for(int i = 0; i < *nrdraws; i++){
    m = draws[i * *ncdraws];
    //freq = draws[i * *ncdraws + 1];
    tau = (int *)my_calloc(m, sizeof(int));
    for(int j = 0; j < m; j++) tau[j] = draws[i * *ncdraws + 2 + j];
    for(int seg = (m-1); seg >= 0; seg--){
      thistau = tau[seg];
      if(seg == 0){
        prevtau = tau[m-1];
      }else{
        prevtau = tau[seg - 1];
      }
      thisStats = (double *)my_calloc(numSumStats, sizeof(double));
      while(thistau != prevtau){
        for(int k = 0; k < numSumStats; k++){
          thisStats[k] += SumStats[k+1][thistau - 1];
        }
        thistau--;
        if(thistau == 0) thistau += *N;
      }

//      Rprintf("%d) m:%d t1:%d t2:%d s:%f n:%f\n", i, m, prevseg, thistau,
//              thisStats[0], thisStats[1]);

      thisSuff = (double *)my_calloc(*nSuffStats, sizeof(double));
      Sufficient_Stats( thisStats, numSumStats, *nSuffStats, Phyp, thisSuff);
      my_free(thisStats);

      for(int k = 0; k < *nSuffStats; k++){
        out[i * *nSuffStats * (*ncdraws - 2) + k * (*ncdraws - 2) + seg] =
          thisSuff[k];
      }

      my_free(thisSuff);
    }
    my_free(tau);
  }

  Free_Summary_Stats(SumStats);

  ESC1:;

  return;
}


//------------------------------------------------------------



void R_init_PeriodCPT_Cfunctions(DllInfo *info){

 R_CMethodDef cMethods[] = {
 {"PeriodCPT_RJMCMC", (DL_FUNC) &PeriodCPT_RJMCMC, 20},
 {"PeriodCPT_EvaluateSufficientStats", (DL_FUNC) &PeriodCPT_EvaluateSufficientStats, 14},
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
// 3: Sufficient Statistic function not identified




