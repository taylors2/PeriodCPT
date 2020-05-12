#ifndef FILE_MCMC
#define FILE_MCMC

#include "PeriodCPT_General.h"
#include "PeriodCPT_MCMCgeneric.h"
#include "PeriodCPT_Cache.h"
#include "PeriodCPT_Proposal.h"


MCMCitem_t * GetInital_tau(int *tau, int *N, int *maxM, int *blank,
                           double *g1, segval_t *g2){
  int m = 0;
  while(tau[m] && m < *maxM) m++;
  return Make_MCMCitem(m, tau, N, g1, g2, 0);
}

void single_step(cache_t **PROPCACHE, int *ncache, int*maxcache,
                  MCMCitem_t *current, double *g1, segval_t *g2,
                  int *minseglen, int *N){

  chain_t *proposals;
  MCMCitem_t *accept;
  double alpha;

  current->j = floor(runif(0,current->m));
  proposals = EvaluateProposal(PROPCACHE, ncache, maxcache, current, minseglen,
                               N, g1, g2);
  accept = sampleProposal(proposals);

  alpha = AcceptanceProb(current, accept, PROPCACHE, ncache, maxcache,
                         g1, g2, minseglen, N);

  if(runif(0,1) < alpha){
    Delete_MCMCitem(current);
    current = accept;
  }else{
    Delete_MCMCitem(accept);
  }
  return;
}

void PushToOutput(int *draw, MCMCitem_t *current, int *blank, int *maxM){
  for(int i = 0; i < current->m; i++){
    draw[i] = current->tau[i];
  }
  for(int i = current->m; i < *maxM; i++){
    draw[i] = *blank;
  }
  return;
}

void PeriodCPT_MCMC(MCMCitem_t *current, int *nsteps, int save,
          int *draw, int *N, int *maxM, int *minseglen, double *g1,
          segval_t *g2, cache_t **PROPCACHE, int *ncache, int *maxcache,
          int *quiet, char **progress_text, int *blank, int *tk, int *err){

  for(int i = 0; i < *nsteps; i++){
    if(*quiet != TRUE) Progress(i, *nsteps, *tk, progress_text);
    single_step(PROPCACHE, ncache, maxcache, current, g1, g2, minseglen, N);
    if(save == TRUE) PushToOutput(&(draw[i * *maxM]), current, blank, maxM);
  }

  return;
}



#endif //FILE_MCMC
