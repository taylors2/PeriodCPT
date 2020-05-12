#ifndef FILE_PROPOSAL
#define FILE_PROPOSAL

#include "PeriodCPT_General.h"
#include "PeriodCPT_MCMCgeneric.h"
#include "PeriodCPT_Cache.h"

int CheckTauOrder(int *tau, int i, int len, int *N){

  if(len == 1){ //M1, no cpts so do nothing
    return i;
  }

  int m, mk;
  mk = *N+1;
  m=0;
  for(int k = 0; k < len; k++){ //check magnitude within [1,N]
    if(tau[k] <= 0){
      tau[k] += *N;
    }else if(tau[k] > *N){
      tau[k] -= *N;
    }
    if(tau[k] < mk){ //find first cpt
      mk = tau[k];
      m = k;
    }
  }

  if(m != 0){ //first cpt is not in first position
    int tautmp[len];
    for(int l = 0;l < len; l++){
      tautmp[l] = tau[(l+m)%len]; //rotate vector
    }
    for(int l = 0; l < len; l++){
      tau[l] = tautmp[l]; //assign back to tau
    }
  }
  i -= m;  //correct position of i
  if(i<0) i+=len;

  return i;
}

void EvaluateProposalsCases(chain_t *proposals, int *tauprop, int *jpos,
                            int mprop, int minseglen, int *N, double *g1,
                            segval_t *g2, int DEATH){
  int st;
  int ipos = *jpos;
  ipos = CheckTauOrder(tauprop, ipos, mprop, N);

  if(DEATH==TRUE){
    st = tauprop[ipos]-minseglen;
  }else if(ipos == 0){
    st = tauprop[mprop-1] + minseglen;
  }else{
    st = tauprop[ipos-1] + minseglen;
  }
  st = ((st-1) % *N) + 1;

  int evaluate = TRUE;
  MCMCitem_t *mcmc;

  while(evaluate == TRUE){
    if(tauprop[ipos]==st) evaluate = FALSE;
    mcmc = Make_MCMCitem(mprop, tauprop, N, g1, g2, ipos);
    Sort_to_chain(proposals, mcmc);
    tauprop[ipos]--;
    if(evaluate == TRUE) ipos = CheckTauOrder(tauprop, ipos, mprop, N);
  }
  *jpos = ipos;
  return;
}

void death_proposals(chain_t *proposals, MCMCitem_t *current,
                     double *g1, segval_t *g2, int *minseglen, int *N){

  //cannot remove any more segments for M1 case
  if(current->m == 1) {
    return;
  }

  if(current->m == 2){ //find pre-defined M1 case in CoreCache.
    int tauprop = 1;
    int ipos = 0;
    //M2->M1 death
    EvaluateProposalsCases(proposals, &tauprop, &ipos, 1, 0, N, g1, g2, TRUE);
    return;
  }

  int *tauprop, ipos, mprop;
  mprop = current->m - 1;
  ipos = current->j;
  tauprop = (int*)calloc(mprop, sizeof(int));
  for(int i = 0; i < mprop; i++) tauprop[i] = current->tau[i + (i>=ipos)];
  if(ipos == mprop) ipos = 0;

  //Move next cpt within minseglen range
  EvaluateProposalsCases(proposals, tauprop, &ipos, mprop,
                         *minseglen-1, N, g1, g2, TRUE);  //general death
  free(tauprop);
  return;
}

void move_proposals(chain_t *proposals, MCMCitem_t *current,
                    double *g1, segval_t *g2, int *minseglen, int *N){

  if(current->m == 1){
    //If M=1, then evaluate case direcly (single scenario)
    int tauprop = 1;
    int ipos = 0;
    //M1 move
    EvaluateProposalsCases(proposals, &tauprop, &ipos, 1, 0, N, g1, g2, TRUE);
    return;
  }

  int *tauprop, mprop, ipos;
  mprop = current->m;
  tauprop = (int*)calloc(mprop, sizeof(int));
  for(int i = 0; i < mprop; i++) tauprop[i] = current->tau[i];
  ipos = current->j;

  if(ipos==(mprop-1)){
    tauprop[ipos] = tauprop[0] - *minseglen;
  }else{
    tauprop[ipos] = tauprop[ipos+1] - *minseglen;
  }
  if(tauprop[ipos]<1) tauprop[ipos] += *N;
  EvaluateProposalsCases(proposals, tauprop, &ipos, mprop,
                         *minseglen, N, g1, g2, FALSE);  //general move
  free(tauprop);

  return;
}


void birth_proposals(chain_t *proposals, MCMCitem_t *current,
                     double *g1, segval_t *g2, int *minseglen, int *N){

  //Get all cases if performing M1->M2 transision
  if(current->m == 1){

    int mprop = current->m + 1;
    int *tauprop, st;
    MCMCitem_t *mcmc, *mcmc2;
    tauprop = (int *)calloc(2, sizeof(int));
    for(tauprop[1] = *minseglen+1; tauprop[1] <= *N; tauprop[1]++){
      st = (tauprop[1] <= (*N - *minseglen)) ? 1 : tauprop[1] + *minseglen - *N ;
      for(tauprop[0] = st; (tauprop[1] - tauprop[0]) >= *minseglen; tauprop[0]++){
        mcmc = Make_MCMCitem(mprop, tauprop, N, g1, g2, 0);
        Sort_to_chain(proposals, mcmc);
        mcmc2 = Copy_MCMCitem(mcmc);
        mcmc2->j = 1;
        Sort_to_chain(proposals, mcmc2);
      }
    }
    free(tauprop);
    return;
  }

  //add to proposal the birth cases
  int mprop = current->m + 1;
  int *tauprop, ipos;
  tauprop = (int*)calloc(mprop,sizeof(int));
  for(int l = 0; l < *minseglen; l++){
    ipos = current->j;
    for(int i = 0; i < (mprop-1); i++){
      tauprop[i+(i>=ipos)] = current->tau[i];
    }
    tauprop[ipos+1] += l;

    //escape if push back too close to next cpt
    if((ipos + 2 == mprop) & ((tauprop[0] + *N - tauprop[mprop-1]) < *minseglen)){
      goto DONE;
    }
    if((ipos + 2 < mprop) & ((tauprop[ipos+2] - tauprop[ipos+1]) < *minseglen)){
      goto DONE;
    }

    tauprop[ipos] = tauprop[ipos+1] - *minseglen; //first position of partition

    //move to next case if insufficient length
    if((ipos == 0) & ((tauprop[0] + *N - tauprop[mprop-1]) < *minseglen)){
      goto HERE;
    }
    if((ipos > 0) & ((tauprop[ipos] - tauprop[ipos - 1]) < *minseglen)){
      goto HERE;
    }

    EvaluateProposalsCases(proposals, tauprop, &ipos, mprop,
                           *minseglen, N, g1, g2, FALSE);  //general birth

    HERE:;
  }
  DONE:;
  free(tauprop);

  return;
}


chain_t *MakeProposals(MCMCitem_t *current, double *g1, segval_t *g2,
                       int *minseglen, int *N){

  chain_t *proposals;
  proposals = Make_Chain();
  death_proposals(proposals, current, g1, g2, minseglen, N);
  move_proposals( proposals, current, g1, g2, minseglen, N);
  birth_proposals(proposals, current, g1, g2, minseglen, N);
  Eval_CMF(proposals);
  return proposals;
}

chain_t *EvaluateProposal(cache_t **PROPCACHE, int *ncache, int*maxcache,
                          MCMCitem_t *current, int *minseglen, int *N,
                          double *g1, segval_t *g2){
  chain_t *proposals;
  proposals = Find_in_Cache(PROPCACHE, ncache, current);
  if(proposals == NULL){
    proposals = MakeProposals(current, g1, g2, minseglen, N);
    Push_To_Cache(PROPCACHE, ncache, maxcache,
                  Make_Cache_Item(Copy_MCMCitem(current), proposals));

  }
  return proposals;
}

MCMCitem_t *sampleProposal(chain_t *proposals){
  double u = runif(0, proposals->last->prob);
  MCMCitem_t *this = proposals->first;
  while(u>this->prob) this = this->next;
  return Copy_MCMCitem(this);
}

double AcceptanceProb(MCMCitem_t *current, MCMCitem_t *accept,
                      cache_t **PROPCACHE, int *ncache, int *cachemax,
                      double *g1, segval_t *g2, int *minseglen, int *N){

  double lZ0, lZ1, lZ2, lmax, alpha;
  chain_t *psetj;

  if(Compare_MCMCitem(accept, current, FALSE)){
    return 1;
  }


  psetj = EvaluateProposal(PROPCACHE, ncache, cachemax, current,
                           minseglen, N, g1, g2);
  lZ0 = log(psetj->last->prob) + psetj->first->value;

  psetj = EvaluateProposal(PROPCACHE, ncache, cachemax, accept,
                           minseglen, N, g1, g2);
  lZ1 = log(psetj->last->prob) + psetj->first->value;

  if(current->m == 1){
    if(accept->j == 1){
      accept->j = 0;
    }else{
      accept->j = 1;
    }
    psetj = EvaluateProposal(PROPCACHE, ncache, cachemax, accept,
                             minseglen, N, g1, g2);
    lZ2 = log(psetj->last->prob) + psetj->first->value;  //Z[kjp]
    lmax = (lZ2>lZ1) ? lZ2 : lZ1;
    alpha = lZ0 - lmax - 2*log(2) + log(exp(lmax - lZ1) + exp(lmax - lZ2));
  }else if(accept->m == 1){
    if(current->j == 1){
      current->j = 0;
    }else{
      current->j = 1;
    }
    psetj = EvaluateProposal(PROPCACHE, ncache, cachemax, accept,
                             minseglen, N, g1, g2);
    lZ2 = log(psetj->last->prob) + psetj->first->value;  //Z[kj]
    lmax = (lZ2>lZ0) ? lZ2 : lZ0;
    alpha = lmax - lZ1 + 2*log(2) - log(exp(lmax - lZ0) + exp(lmax - lZ2));
  }else{
    alpha = lZ0 - lZ1 + log(current->m) - log(accept->m);
  }

  alpha = (alpha>0) ? 1 : exp(alpha);

  return alpha;
}


#endif //FILE_PROPOSAL
