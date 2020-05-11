#include "PeriodCPT_General.h"
#include "PeriodCPT_LOOKUP.h"
//#include "PeriodCPT_cache.h"

MCMCitem_t* Make_MCMCitem(int m, int *tau, int *N, double *g1, segval_t *g2, int j){
  MCMCitem_t *mcmc = calloc(1, sizeof(MCMCitem_t));
  mcmc->tau = (int *)calloc(m, sizeof(MCMCitem_t));
  for(int i = 0; i < m; i++) mcmc->tau[i] = tau[i];
  mcmc->j = j;
  mcmc->m = m;
  mcmc->value = Eval_at_CPT(tau, m, N, g1, g2);
  mcmc->prob = 0;
  mcmc->prev = NULL;
  mcmc->next = NULL;
  return mcmc;
}

void Delete_MCMCitem(MCMCitem_t *mcmc){
  if(mcmc->prev != 0) mcmc->prev = mcmc->next;
  if(mcmc->next != 0) mcmc->next = mcmc->prev;
  free(mcmc->tau);
  free(mcmc);
  return;
}

chain_t * Make_Chain(){
  chain_t *chain;
  chain = (chain_t*)calloc(1, sizeof(chain_t));
  chain->first  = NULL;
  chain->last   = NULL;
  chain->length = 0;
  return chain;
}

void Delete_Chain(chain_t *chain){
  MCMCitem_t *mcmc;
  while(chain->first != NULL){
    mcmc = chain->first;
    chain->first = mcmc->next;
    Delete_MCMCitem(mcmc);
    chain->length--;
  }
  free(chain);
  return;
}



MCMCitem_t * GetInital_tau(int *tau, int *N, int *maxM, int *blank,
                           double *g1, segval_t *g2){
  int m = 0;
  while(tau[m] && m < *maxM) m++;
  return Make_MCMCitem(m, tau, N, g1, g2, 0);
}


//--------------------------------------------------------------------------

cache_t *Make_Cache_Item(MCMCitem_t *generator, chain_t *chain){
  cache_t *p   = (cache_t *)calloc(1, sizeof(cache_t));
  p->generator = generator;
  p->chain     = chain;
  return p;
}

void Delete_Cache_Item(cache_t *item){
  if(item != 0){
    Delete_MCMCitem(item->generator);
    Delete_Chain(item->chain);
  }
  return;
}

cache_t **Make_Cache_List(int maxlen){
  cache_t **p = (cache_t **)calloc(maxlen,sizeof(cache_t *));
  return p;
}

void Delete_Cache_List(cache_t **list, int maxlen){
  for(int i = 0; i < maxlen; i++){
    Delete_Cache_Item(list[i]);
    list[i] = 0;
  }
  free(list);
}

int Compare_Cache_Item(cache_t *cache, MCMCitem_t *mcmc){
  if(cache->generator->m != mcmc->m) return 0;
  if(cache->generator->j != mcmc->j) return 0;
  for(int j = 0; j < mcmc->m; j++){
    if(cache->generator->tau[j] != mcmc->tau[j]) return 0;
  }
  return 1;
}

cache_t *Find_in_Cache(cache_t **cache, int *n, MCMCitem_t *mcmc){
  cache_t *match;
  match = NULL;
  for(int i = 0; (i < *n) && (match == NULL); i++){
    if(Compare_Cache_Item(cache[i], mcmc)){
      match = cache[i];
    }
  }
  return match;
}

//--------------------------------------------------------------------------
chain_t *EvaluateProposal(cache_t **PROPCACHE, int *ncache, int*maxcache,
                          MCMCitem_t *current, int *minseglen, int *N,
                          double *g1, segval_t *g2){
  chain_t *proposals;
  proposals = Find_in_Cache(PROPCACHE, ncache, current);
  if(proposals == NULL){
    proposals = MakeProposals(current, g1, g2, minseglen, N);
    PushToList(PROPCACHE, current, proposals);
  }
  return proposals;
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
    free(current->tau);
    current->m = accept->m;
    current->value = accept->value;
    current->j = accept->j;
    current->tau = (int*)calloc(current->m,sizeof(int));
    int i;
    for(int i = 0; i < current->m; i++) current->tau[i] = accept->tau[i];
  }
  DeleteMCMCitem(accept);

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
          int *quiet, char **progress_text, int *blank, int *err){

  for(int i = 0; i < *nsteps; i++){
    if(*quiet != TRUE) Progress(i, *nsteps, progress_text);
    single_step(PROPCACHE, ncache, maxcache, current, g1, g2, minseglen, N);
    if(save == TRUE) PushToOutput(&(draw[i * *maxM]), current, blank, maxM);
  }

  return;
}


