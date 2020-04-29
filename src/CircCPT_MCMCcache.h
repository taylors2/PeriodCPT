/*

#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\stdio.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\stdlib.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\math.h>
#include <C:\Program Files\R\R-3.4.3\include\R.h>
#include <C:\Program Files\R\R-3.4.3\include\Rmath.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\time.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\stdint.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\string.h>
#include"CircCPT_General.c"
#include"CircCPT_LOOKUP.c"

 */



//--allocate memory

chain_t *MakeMCMCchain(){
int fnname = 200;
  profile(TRUE, fnname);

  chain_t *CHAIN;
  CHAIN = (chain_t *)my_calloc(1, sizeof(chain_t));
  CHAIN->first = NULL;
  CHAIN->last = NULL;
  CHAIN->length = 0;

  profile(FALSE, fnname);
  return CHAIN;
}

MCMCitem_t* MakeMCMCitem(int m, int *tau, int *N, double *g1, struct G2 *g2, int j){
  int fnname = 201;
  profile(TRUE, fnname);

  MCMCitem_t *mcmc;
  int i;
  mcmc = (MCMCitem_t *)my_calloc(1, sizeof(MCMCitem_t));
  mcmc->tau = (int *)my_calloc(m, sizeof(MCMCitem_t));
  for(i=0;i<m;i++) mcmc->tau[i] = tau[i];
  mcmc->j = j;
  mcmc->m = m;
  mcmc->value = Eval_at_CPT(tau, m, N, g1, g2);
  mcmc->prob = 0;
  mcmc->next = NULL;
  mcmc->prev = NULL;

  profile(FALSE, fnname);
  return mcmc;
}

MCMCitem_t *CopyMCMCitem(MCMCitem_t *mcmc){
  int fnname = 202;
  profile(TRUE, fnname);

  MCMCitem_t *new;
  int i;
  new = (MCMCitem_t *)my_calloc(1, sizeof(MCMCitem_t));
  new->tau = (int *)my_calloc(mcmc->m, sizeof(MCMCitem_t));
  for(i=0; i<mcmc->m; i++) new->tau[i] = mcmc->tau[i];
  new->j = mcmc->j;
  new->m = mcmc->m;
  new->value = mcmc->value;
  new->prob = 0;
  new->next = NULL;
  new->prev = NULL;

  profile(FALSE, fnname);
  return new;
}

//--free memeory

void DeleteMCMCitem(MCMCitem_t *mcmc){
  int fnname = 203;
  profile(TRUE, fnname);

  my_free(mcmc->tau);
  my_free(mcmc);

  profile(FALSE, fnname);
  return;
}

void DeleteMCMCchain(chain_t *CHAIN){
  int fnname = 204;
  profile(TRUE, fnname);

  MCMCitem_t *this, *next;
  this = CHAIN->first;
  while(this != NULL){
    next = this->next;
    DeleteMCMCitem(this);
    this = next;
  }
  my_free(CHAIN);

  profile(FALSE, fnname);
  return;
}

//-- Push & Pop

void PushToChain(chain_t *CHAIN, MCMCitem_t *mcmc){ //append mcmc to CHAIN
  int fnname = 205;
  profile(TRUE, fnname);

  if(CHAIN->first == NULL){
    CHAIN->first = mcmc;
    CHAIN->last = mcmc;
    mcmc->prev = NULL;
    mcmc->next = NULL;
  }else{
    CHAIN->last->next = mcmc;
    mcmc->prev = CHAIN->last;
    CHAIN->last = mcmc;
    mcmc->next = NULL;
  }
  CHAIN->length++;

  profile(FALSE, fnname);
  return;
}

void BubbleToChain(chain_t *CHAIN, MCMCitem_t *this){ //append mcmc to CHAIN
  int fnname = 206;
  profile(TRUE, fnname);

  if(CHAIN->first == NULL){
    CHAIN->first = this;
    CHAIN->last = this;
    CHAIN->length = 1;
    profile(FALSE, fnname);
    return;
  }

  if(CHAIN->first->value > this->value + 20){
    //prob of this will be ~0, push straight to end.
    this->prev = CHAIN->last;
    this->next = NULL;
    this->prev->next = this;
    CHAIN->last = this;
    CHAIN->length++;
    profile(FALSE, fnname);
    return;
  }

  MCMCitem_t *place;
  int point = FALSE;
  place = CHAIN->first;
  while(point == FALSE){
    if(place==NULL){
      point = TRUE;
    }else if(place->value <= this->value){
      point = TRUE;
    }else{
      place = place->next;
    }
  }
  if(place == NULL){ //push to end
    this->prev = CHAIN->last;
    this->next = NULL;
    this->prev->next = this;
    CHAIN->last = this;
  }else if(place->prev == NULL){ //push to front
    place->prev = this;
    this->next = place;
    this->prev = NULL;
    CHAIN->first = this;
  }else{
    place->prev->next = this;
    this->prev = place->prev;
    place->prev = this;
    this->next = place;
  }
  CHAIN->length++;

  profile(FALSE, fnname);
  return;
}

void PopFromChain(chain_t *CHAIN, MCMCitem_t *mcmc, int FREE){
  int fnname = 207;
  profile(TRUE, fnname);

  if(CHAIN->last != mcmc){
    if(CHAIN->first != mcmc){  //somewhere in the middle
      mcmc->prev->next = mcmc->next;
      mcmc->next->prev = mcmc->prev;
    }else{  //is the first in chain
      CHAIN->first = mcmc->next;
      CHAIN->first->prev = NULL;
    }
  }else{
    if(CHAIN->first != mcmc){ //is the last in chain
      CHAIN->last = mcmc->prev;
      CHAIN->last->next = NULL;
    }else{  //is the only item in chain
      CHAIN->first = NULL;
      CHAIN->last = NULL;
    }
  }
  mcmc->prev = NULL;
  mcmc->next = NULL;
  CHAIN->length--;

  if(FREE == TRUE) DeleteMCMCitem(mcmc);
  profile(FALSE, fnname);
  return;
}

//-- Clear chain to retain last case

void ClearChain(chain_t *CHAIN){
  int fnname = 208;
  profile(TRUE, fnname);

  if(CHAIN->first == CHAIN->last){
    profile(FALSE, fnname);
    return;
  }

  MCMCitem_t *this, *next;
  CHAIN->last->prev->next = NULL;
  CHAIN->last->prev = NULL;
  this = CHAIN->first;
  while(this!=NULL){
    next = this->next;
    DeleteMCMCitem(this);
    this = next;
  }
  CHAIN->first = CHAIN->last;
  CHAIN->length = 1;

  profile(FALSE, fnname);
  return;
}

//-- Print

void PrintMCMCchain(chain_t *CHAIN, int index){
  int fnname = 209;
  profile(TRUE, fnname);

  MCMCitem_t *mcmc;
  mcmc = CHAIN->first;
  int maxm;
  maxm = 0;
  while(mcmc != NULL){
    if(maxm < mcmc->m) maxm = mcmc->m;
    mcmc = mcmc->next;
  }

  FILE *f;
  int j;
  char fname[100];
  sprintf(fname,"%d_MCMClist.txt", index);

  f = fopen(fname, "w");
  fprintf(f,"me, prev, next, id, m, j, value, prob");
  for(j=0; j<maxm; j++) fprintf(f,", tau%d", j+1);
  fprintf(f,"\n");

  mcmc = CHAIN->first;
  int i = 0;
  while(i < CHAIN->length){
    fprintf(f,"%p, ", (void *)mcmc);
    if(mcmc->prev == NULL){
      fprintf(f,"-1, ");
    }else{
      fprintf(f,"%p, ", (void *)mcmc->prev);
    }
    if(mcmc->next == NULL){
      fprintf(f,"-1, ");
    }else{
      fprintf(f,"%p, ", (void *)mcmc->next);
    }
    fprintf(f,"%d, %d, %d, %f, %f", i, mcmc->m, mcmc->j, mcmc->value,mcmc->prob);
    for(j=0;j<maxm;j++){
      if(j<mcmc->m){
        fprintf(f,", %d", mcmc->tau[j]);
      }else{
        fprintf(f,", -1");
      }
    }
    fprintf(f,"\n");
    mcmc = mcmc->next;
    i++;
  }
  fclose(f);

  profile(FALSE, fnname);
  return;
}

// -- Export

int ReturnChain(chain_t *CHAIN, int *draw, int offset, int *ndraw, int *err, double id){
  int fnname = 210;
  profile(TRUE, fnname);

  int i, j;
  MCMCitem_t *this;
  this = CHAIN->first;
  i = offset;
  while(this != NULL){
    if(*err != 0) break;

    //write cpt vector
    for(j=0; j < this->m; j++){
      if(i<*ndraw){
        draw[i++] = this->tau[j];
      }else{
        *err = 2;
        break;
      }
    }

    //push -id at end of vector
    if(i<*ndraw){
      draw[i++] = -id;
    }else{
      *err = 2;
    }

    this = this->next;
  }

  profile(FALSE, fnname);
  return i;
}

void ExportDraws(chain_t *chain1, chain_t *chain2, int *draw, int *ndraw, int *err){
  int fnname = 211;
  profile(TRUE, fnname);
  int i=0;
  i = ReturnChain(chain1, draw, i, ndraw, err, 1.0);
  i = ReturnChain(chain2, draw, i, ndraw, err, 2.0);
  profile(FALSE, fnname);
  return;
}

void EvalProbs(chain_t *chain){
  int fnname = 212;
  profile(TRUE, fnname);
  double max, sum;
  MCMCitem_t *this;
  max = chain->first->value;

  //This can be commented if chain is sorted as chain->first is maximum
  //this = chain->first;
  //while(this!=NULL){
  //  if(max < this->value) max = this->value;
  //  this = this->next;
  //}

  this = chain->first;
  sum = 0.0;
  while(this!=NULL){
    sum += exp(this->value - max);
    this->prob = sum;
    this = this->next;
  }

  profile(FALSE, fnname);
  return;
}


int CompareMCMCitem(MCMCitem_t *test, MCMCitem_t *instance){
  int fnname = 213;
  profile(TRUE, fnname);
  int i;
  if(test->m != instance->m) goto DIFF;
  if(test->j != instance->j) goto DIFF;
  for(i=0;i<test->m;i++){
    if(test->tau[i] != instance->tau[i]) goto DIFF;
  }
  profile(FALSE, fnname);
  return TRUE;

  DIFF:;
  profile(FALSE, fnname);
  return TRUE;  //ERROR?!
}
