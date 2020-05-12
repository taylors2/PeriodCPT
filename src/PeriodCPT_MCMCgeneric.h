#ifndef FILE_MCMCGENERIC
#define FILE_MCMCGENERIC

#include "PeriodCPT_General.h"
#include "PeriodCPT_LOOKUP.h"

MCMCitem_t* Make_MCMCitem(int m, int *tau, int *N,
                          double *g1, segval_t *g2, int j){
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


void Push_to_chain(chain_t *chain, MCMCitem_t *item){
  chain->last->next = item;
  item->next = NULL;
  item->prev = chain->last->next;
  chain->last = item;
  chain->length++;
}

void Sort_to_chain(chain_t *chain, MCMCitem_t *item){
  MCMCitem_t *check;
  check = chain->first;
  int located = FALSE;
  while(located == FALSE && check != NULL){
    if(check->value > item->value){
      check = check->next;
    }else{
      located = TRUE;
    }
  }
  if(located == FALSE){   //item has the smallest value
    Push_to_chain(chain, item);
  }else if(check == chain->first){  //item has the largest value
    item->prev = NULL;
    item->next = check;
    check->prev = item;
    chain->first = item;
    chain->length++;
  }else{  //item appears somewhere in the middle of chain
    item->prev = check->prev;
    item->next = check;
    check->prev->next = item;
    check->prev = item;
    chain->length++;
  }
  return;
}

void Pop_from_chain(chain_t *chain, MCMCitem_t *item){
  if(item->prev != NULL) item->prev->next = item->next;
  if(item->next != NULL) item->next->prev = item->prev;
  if(chain->first == item) chain->first = item->next;
  if(chain->last == item) chain->last = item->prev;
  chain->length--;
  item->prev = NULL;
  item->next = NULL;
  return;
}

MCMCitem_t *Copy_MCMCitem(MCMCitem_t *mcmc){
  MCMCitem_t *newmcmc;
  newmcmc = (MCMCitem_t *)calloc(1, sizeof(MCMCitem_t));
  newmcmc->j        = mcmc->j;
  newmcmc->m        = mcmc->m;
  newmcmc->next     = NULL;
  newmcmc->prev     = NULL;
  newmcmc->prob     = mcmc->prob;
  newmcmc->value    = mcmc->value;
  newmcmc->tau      = (int *)calloc(mcmc->m, sizeof(int));
  for(int i = 0; i < mcmc->m; i++){
    newmcmc->tau[i] = mcmc->tau[i];
  }
  return newmcmc;
}

void Delete_Chain(chain_t *chain){
  MCMCitem_t *mcmc;
  while(chain->first != NULL){
    mcmc = chain->first;
    Pop_from_chain(chain, mcmc);
    Delete_MCMCitem(mcmc);
  }
  free(chain);
  return;
}


void Eval_CMF(chain_t *chain){
  //Assume chain is sorted
  double max, sum;
  MCMCitem_t *this;
  max = chain->first->value;

  this = chain->first;
  sum = 0.0;
  while(this!=NULL){
    sum += exp(this->value - max);
    this->prob = sum;
    this = this->next;
  }
  return;
}

int Compare_MCMCitem(MCMCitem_t *test, MCMCitem_t *instance, int checkj){
  if(test->m != instance->m) return 0;
  if((checkj == TRUE) && (test->j != instance->j)) return 0;
  for(int j = 0; j < instance->m; j++){
    if(test->tau[j] != instance->tau[j]) return 0;
  }
  return 1;
}

#endif //FILE_MCMCGENERIC
