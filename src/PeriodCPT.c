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
#include "PeriodCPT_modepcpt.h"

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
    int     *draw,       //Vectorised MCMC output (append each chain)
    int     *Eval_Fits,  //Logical, should mode and fit estimates be calculated?
    int     *mode_pcpt,  //pcpt mode output
    double  *mode_params,//Segment parameter posterior mode output
    int     *nsegparams, //Number of parameters per segment
    double  *fits,       //Fit metrics output
    int     *nfits       //Number of fit metrics
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

  if(*Eval_Fits == 1){
    int len = *nchains * *niter;
    Mode_pcpt(draw, maxM, &len, blank, N, g1, g2, mode_pcpt);
    Mode_Fit_Calcuation(data, time, n, N, Pdist, Phyp, mode_pcpt, blank, maxM, err,
                        fits, nfits, mode_params, nsegparams);
  }

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
    char   **Pdist,      //Name of sampling distribution
    double  *Phyp,       //Hyper-parameter(s) for segment parameter priors
    int     *draws,      //m freq tau1, ...
    int     *ncdraws,    //length of draw per case
    int     *nrdraws,    //number of cases in draw
    int     *nSuffStats, //number of (posterior) sufficient statistics
    int     *blank,      //Holding number of blank or NA spaces in inits and draw
    int     *err,        //Error flag
    double  *out         //Vectorised output
){

  *err = 0;

  //Grab distribution related functions
  Samp_Dist_Ptr *Samp_Dist = NULL;
  Summary_Stats_Ptr *Summary_Stats = NULL;
  Sufficient_Stats_Ptr *Sufficient_Stats = NULL;
  Get_Pprior(Pdist, &Samp_Dist, err);
  Get_SumStats_FN(Pdist, &Summary_Stats, err);
  Get_SuffStats_Function(Pdist, &Sufficient_Stats, err);
  if(*err != 0){
    goto ESC1;
  }
  //flood output with blank character
  for(int i = 0; i < *nSuffStats * (*ncdraws) * (*nrdraws); i++){
    out[i] = *blank;
  }

  double **SumStats = Summary_Stats(data, time, n, N);
  int numSumStats = (int)SumStats[0][0];

  int thistau, prevtau;
  double *thisStats, *thisSuff;
  MCMCitem_t *current;

  for(int i = 0; i < *nrdraws; i++){
    int m = 0;
    for(int j=0; j<*ncdraws; j++){
      if(draws[i * *ncdraws + j] != *blank) m++;
    }
    current = Make_blank_MCMCitem(m, &(draws[i * *ncdraws]));

    for(int seg = (current->m-1); seg >= 0; seg--){
      thistau = current->tau[seg];
      if(seg == 0){
        prevtau = current->tau[current->m - 1];
      }else{
        prevtau = current->tau[seg - 1];
      }

      thisStats = (double *)my_calloc(numSumStats, sizeof(double));

      while(thistau != prevtau){
        for(int k = 0; k < numSumStats; k++){
          thisStats[k] += SumStats[k+1][thistau - 1];
        }
        thistau--;
        if(thistau == 0) thistau += *N;
      }

      thisSuff = (double *)my_calloc(*nSuffStats, sizeof(double));
      Sufficient_Stats( thisStats, numSumStats, *nSuffStats, Phyp, thisSuff);
      my_free(thisStats);

      for(int k = 0; k < *nSuffStats; k++){
        out[i * *nSuffStats * (*ncdraws) + k * (*ncdraws) + seg] =
          thisSuff[k];
      }

      my_free(thisSuff);
    }
    Delete_MCMCitem(current);
  }

  Free_Summary_Stats(SumStats);

  ESC1:;

  return;
}

//-----------------------------------------------------------

extern void Tally_pcpt(int *chains, int *nrow, int *ncol, int *blank, int *tally){
  //Create table of unique pcpt samples, using the prob slot at counter
  chain_t *unique = Make_Chain();
  MCMCitem_t *current = NULL;
  List_unique_samples(chains, ncol, nrow, blank, unique);

  current = unique->first;
  while(current != NULL){
    tally[current->j] = (int)current->prob;
    current = current->next;
  }
  Delete_Chain(unique);
  return;
}



//------------------------------------------------------------



void R_init_PeriodCPT_Cfunctions(DllInfo *info){

 R_CMethodDef cMethods[] = {
 {"PeriodCPT_RJMCMC", (DL_FUNC) &PeriodCPT_RJMCMC, 26},
 {"PeriodCPT_EvaluateSufficientStats", (DL_FUNC) &PeriodCPT_EvaluateSufficientStats, 13},
 {"Tally_pcpt", (DL_FUNC) &Tally_pcpt, 5},
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
// 3: Summary Statistic function not identified
// 4: Sufficient Statistic function not identified
// 5: Segment parameter mode function not identified
// 6: Fit metric calculation function not identified
//WARNINGS
// 101: Segment param mode is not unique. Consider changing prior hyper-params.
// 102: Small sample size on segment for reliable fit metric. Consider increasing minimum segment length.





