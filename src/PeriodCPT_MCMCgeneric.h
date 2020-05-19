#ifndef FILE_MCMCGENERIC
#define FILE_MCMCGENERIC

#include "PeriodCPT_General.h"
#include "PeriodCPT_LOOKUP.h"

MCMCitem_t* Make_MCMCitem(int m, int *tau, int *N,
                          double *g1, segval_t *g2, int j){
  //Rprintf("Make_MCMCitem\n");
  MCMCitem_t *mcmc = my_calloc(1, sizeof(MCMCitem_t));
  mcmc->tau = (int *)my_calloc(m, sizeof(MCMCitem_t));
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
  //Rprintf("Delete_MCMCitem\n");
  my_free(mcmc->tau);
  my_free(mcmc);
  return;
}

chain_t * Make_Chain(){
  //Rprintf("Make_Chain\n");
  chain_t *chain;
  chain = (chain_t*)my_calloc(1, sizeof(chain_t));
  chain->first  = NULL;
  chain->last   = NULL;
  chain->length = 0;
  return chain;
}


void Push_to_chain(chain_t *chain, MCMCitem_t *item){
  //Rprintf("Push_to_chain\n");
  if(chain->first == NULL){
    chain->first = item;
    chain->last = item;
    item->next = NULL;
    item->prev = NULL;
    chain->length++;
  }else{
    chain->last->next = item;
    item->next = NULL;
    item->prev = chain->last;
    chain->last = item;
    chain->length++;
  }
  return;
}

void Chain_Insert_Before_Case(MCMCitem_t *item, MCMCitem_t *this, chain_t *chain){
  //Rprintf("Chain_Insert_Before_Case\n");
  //Push 'item' before 'this'
  if(this == NULL){ //this is undefined, so push item to end
    Push_to_chain(chain, item);
    return;
  }
  if(this->prev != NULL) this->prev->next = item;
  item->prev = this->prev;
  item->next = this;
  this->prev = item;
  if(item->prev == NULL) chain->first = item;
  if(item->next == NULL) chain->last = item;  //item not at end, soshould not be called!
  chain->length++;
  return;
}


void Sort_to_chain(chain_t *chain, MCMCitem_t *item){
  //Rprintf("Sort_to_chain\n");
  MCMCitem_t *this;
  this = chain->first;
  int search = TRUE;
  while(search == TRUE){
    if(this == NULL){ //Hit chain end (item has smallest value)
      Push_to_chain(chain, item);
      return;
    }else if(this->value < item->value){
      search = FALSE;
    }else{
      this = this->next;
    }
  }
  Chain_Insert_Before_Case(item, this, chain);
  return;
}

void Pop_from_chain(chain_t *chain, MCMCitem_t *item){
  //Rprintf("Pop_from_chain\n");
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
  //Rprintf("Copy_MCMCitem\n");
  MCMCitem_t *newmcmc;
  newmcmc = (MCMCitem_t *)my_calloc(1, sizeof(MCMCitem_t));
  newmcmc->j        = mcmc->j;
  newmcmc->m        = mcmc->m;
  newmcmc->next     = NULL;
  newmcmc->prev     = NULL;
  newmcmc->prob     = mcmc->prob;
  newmcmc->value    = mcmc->value;
  newmcmc->tau      = (int *)my_calloc(mcmc->m, sizeof(int));
  for(int i = 0; i < mcmc->m; i++){
    newmcmc->tau[i] = mcmc->tau[i];
  }
  return newmcmc;
}

void Delete_Chain(chain_t *chain){
  //Rprintf("Delete_Chain\n");
  MCMCitem_t *mcmc;
  while(chain->first != NULL){
    mcmc = chain->first;
    Pop_from_chain(chain, mcmc);
    Delete_MCMCitem(mcmc);
  }
  my_free(chain);
  return;
}

void Eval_CMF(chain_t *chain){
  //Rprintf("Eval_CMF\n");
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
  //Rprintf("Compare_MCMCitem\n");
  if(test->m != instance->m) return 0;
  if((checkj == TRUE) && (test->j != instance->j)) return 0;
  for(int j = 0; j < instance->m; j++){
    if(test->tau[j] != instance->tau[j]) return 0;
  }
  return 1;
}

void Swap_MCMCitems(MCMCitem_t *a, MCMCitem_t *b){
  //Rprintf("Swap_MCMCitems\n");
  MCMCitem_t tmp = *a;
  *a = *b;
  *b = tmp;
  return;
}

MCMCitem_t * Find_in_Chain(chain_t *chain, MCMCitem_t *mcmc){
  if(chain->first == NULL) return NULL;
  MCMCitem_t *this = chain->first;
  while(this != NULL){
    if(Compare_MCMCitem(mcmc, this, FALSE) == 1){
      return this;
    }
  }
  return NULL;
}

#endif //FILE_MCMCGENERIC
