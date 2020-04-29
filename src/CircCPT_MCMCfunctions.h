/* 

#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\stdio.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\stdlib.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\math.h>
#include <C:\Program Files\R\R-3.4.3\include\R.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\time.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\stdint.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\string.h>
#include"CircCPT_General.c"
#include"CircCPT_MCMCcache.c"
#include"CircCPT_PROPfunctions.c"

 */


int CheckTauOrder2(int *tau, int i, int len, int *N){
  int fnname = 300;
  profile(TRUE, fnname);
  
  if(len == 1){ //M1, no cpts so do nothing
    profile(FALSE, fnname);
    return i;
  }
  
  int m, mk, k;
  mk = *N+1;
  m=0;
  for(k=0;k<len;k++){ //check magnitude within [1,N]
    if(tau[k]<=0){
      tau[k] += *N;
    }else if(tau[k]>*N){
      tau[k] -= *N;
    }
    if(tau[k] < mk){ //find first cpt
      mk = tau[k];
      m = k;
    }
  }
  
  if(m!=0){ //first cpt is not in first position
    int l, tautmp[len];
    for(l=0;l<len;l++){
      tautmp[l] = tau[(l+m)%len]; //rotate vector
    }
    for(l=0;l<len;l++){
      tau[l] = tautmp[l]; //assign back to tau
    }
  }
  i -= m;  //correct position of i
  if(i<0) i+=len;

  profile(FALSE, fnname);
  return i;
}


void EvaluateProposalsCases(chain_t *proposals, int *tauprop, int *jpos, int mprop,
                            int minseglen, int *N, double *g1, struct G2 *g2, int DEATH){
  int fnname = 301;
  profile(TRUE, fnname);

  int st;
  int ipos = *jpos;
  ipos = CheckTauOrder2(tauprop, ipos, mprop, N);

  if(DEATH==TRUE){
    st = tauprop[ipos]-minseglen;
  }else if(ipos == 0){
    st = tauprop[mprop-1] + minseglen;
  }else{
    st = tauprop[ipos-1] + minseglen;
  }

  if(st>*N){
    st -= *N;
  }else if(st < 1){
    st += *N;
  }

  int evaluate = TRUE;
  MCMCitem_t *mcmc;

  while(evaluate == TRUE){
    if(tauprop[ipos]==st) evaluate = FALSE;
    mcmc = MakeMCMCitem(mprop, tauprop, N, g1, g2, ipos);
    BubbleToChain(proposals, mcmc);
    tauprop[ipos]--;
    if(evaluate == TRUE) ipos = CheckTauOrder2(tauprop, ipos, mprop, N);
  }
  *jpos = ipos;
  profile(FALSE, fnname);
  return;
}

void death_proposals(chain_t *proposals, MCMCitem_t *current, double *g1, struct G2 *g2, int *minseglen, int *N){

  int fnname = 302;
  profile(TRUE, fnname);
  //cannot remove any more segments for M1 case
  if(current->m == 1) {
    profile(FALSE, fnname);
    return;  
  }

  if(current->m ==2){ //find pre-defined M1 case in CoreCache.
    int tauprop = 1;
    int ipos = 0;
    EvaluateProposalsCases(proposals, &tauprop, &ipos, 1, 0, N, g1, g2, TRUE);  //M2->M1 death
    profile(FALSE, fnname);
    return;  
  }

  int *tauprop, i, ipos, mprop;
  mprop = current->m - 1;
  ipos = current->j;
  tauprop = (int*)my_calloc(mprop, sizeof(int));
  for(i=0; i < mprop; i++) tauprop[i] = current->tau[i + (i>=ipos)];
  if(ipos == mprop) ipos = 0;

  //Move next cpt within minseglen range
  EvaluateProposalsCases(proposals, tauprop, &ipos, mprop, *minseglen-1, N, g1, g2, TRUE);  //general death
  my_free(tauprop);
  profile(FALSE, fnname);
  return;
}
  
void move_proposals(chain_t *proposals, MCMCitem_t *current, double *g1, struct G2 *g2, int *minseglen, int *N){

  int fnname = 303;
  profile(TRUE, fnname);

  if(current->m == 1){
    //If M=1, then evaluate case direcly (single scenario)
    int tauprop = 1;  
    int ipos = 0;
    EvaluateProposalsCases(proposals, &tauprop, &ipos, 1, 0, N, g1, g2, TRUE);  //M1 move
    profile(FALSE, fnname);
    return;  
  }

  int *tauprop, mprop, i, ipos;
  mprop = current->m;
  tauprop = (int*)my_calloc(mprop, sizeof(int));
  for(i=0;i<mprop;i++) tauprop[i] = current->tau[i];
  ipos = current->j;
  
  if(ipos==(mprop-1)){
    tauprop[ipos] = tauprop[0] - *minseglen;
  }else{
    tauprop[ipos] = tauprop[ipos+1] - *minseglen;
  }
  if(tauprop[ipos]<1) tauprop[ipos]+=*N;
//  ipos = CheckTauOrder2(tauprop, ipos, mprop, N);
  EvaluateProposalsCases(proposals, tauprop, &ipos, mprop, *minseglen, N, g1, g2, FALSE);  //general move
  my_free(tauprop);

  profile(FALSE, fnname);
  return;
}


void birth_proposals2(chain_t *proposals, MCMCitem_t *current, double *g1, struct G2 *g2, int *minseglen, int *N){
  int fnname = 304;
  profile(TRUE, fnname);
  
  //Get all cases if performing M1->M2 transision
  if(current->m == 1){
    
    int mprop = current->m + 1;
    int *tauprop, st;
    MCMCitem_t *mcmc, *mcmc2;
    tauprop = (int *)my_calloc(2,sizeof(int));
    for(tauprop[1]=*minseglen+1; tauprop[1]<=*N; tauprop[1]++){
      st = (tauprop[1] <= (*N-*minseglen)) ? 1 : tauprop[1]+*minseglen-*N ;
      for(tauprop[0]=st; (tauprop[1]-tauprop[0])>= *minseglen; tauprop[0]++){
        mcmc = MakeMCMCitem(mprop, tauprop, N, g1, g2, 0);
        BubbleToChain(proposals, mcmc);
        mcmc2 = CopyMCMCitem(mcmc);
        mcmc2->j = 1;
        BubbleToChain(proposals, mcmc2);
      }
    }
    my_free(tauprop);
    
    profile(FALSE, fnname);
    return;
  }
  
  //add to proposal the birth cases
  int mprop = current->m + 1;
  int *tauprop, i, ipos, l;
  tauprop = (int*)my_calloc(mprop,sizeof(int));
  for(l=0; l<*minseglen; l++){
    ipos = current->j;
    for(i=0;i<(mprop-1);i++) tauprop[i+(i>=ipos)] = current->tau[i];
    tauprop[ipos+1] += l;

    //escape if push back too close to next cpt
    if((ipos+2 == mprop) & ((tauprop[0]+*N-tauprop[mprop-1])<*minseglen)) goto DONE;
    if((ipos+2 < mprop) & ((tauprop[ipos+2]-tauprop[ipos+1])<*minseglen)) goto DONE;

    tauprop[ipos] = tauprop[ipos+1] - *minseglen; //first position of partition

    //move to next case if insufficient length
    if((ipos==0) & ((tauprop[0]+*N-tauprop[mprop-1]) < *minseglen)) goto HERE;
    if((ipos>0) & ((tauprop[ipos]-tauprop[ipos-1]) < *minseglen)) goto HERE;

    EvaluateProposalsCases(proposals, tauprop, &ipos, mprop, *minseglen, N, g1, g2, FALSE);  //general birth
    HERE:;
  }
  DONE:;
  my_free(tauprop);
  
  profile(FALSE, fnname);
  return; 
}

//-- RJMCMC steps

chain_t *MakeProposals(MCMCitem_t *current, double *g1, struct G2 *g2, int *minseglen, int *N){

  int fnname = 305;
  profile(TRUE, fnname);

  chain_t *proposals;
  proposals = MakeMCMCchain();
  death_proposals(proposals, current, g1, g2, minseglen, N);
  move_proposals( proposals, current, g1, g2, minseglen, N);
  birth_proposals2(proposals, current, g1, g2, minseglen, N);
  EvalProbs(proposals);
  profile(FALSE, fnname);
  return proposals;
}


chain_t *EvaluateProposal(list_t *PROPCACHE, MCMCitem_t *current, 
                          int *minseglen, int *N, double *g1, struct G2 *g2){
  int fnname = 314;
  profile(TRUE, fnname);
  
  chain_t *proposals;
  proposals = FindCase(PROPCACHE, current);
  if(proposals == NULL){
    proposals = MakeProposals(current, g1, g2, minseglen, N);
    PushToList(PROPCACHE, current, proposals);
  }

  profile(FALSE, fnname);
  return proposals;
}

MCMCitem_t *sampleProposal(chain_t *proposals){
  int fnname = 306;
  profile(TRUE, fnname);
  double u = runif(0, proposals->last->prob);
  MCMCitem_t *this, *out;
  this = proposals->first;
  while(u>this->prob){
    this = this->next;
  }
  out = CopyMCMCitem(this);

  profile(FALSE, fnname);
  return out;
}


double AcceptanceProb(MCMCitem_t *current, MCMCitem_t *accept,
    list_t *PROPCACHE, double *g1, struct G2 *g2, int *minseglen, int *N){
  int fnname = 314;
  profile(TRUE, fnname);

  int k, SAME;
  double lZ0, lZ1, lZ2, lmax, alpha;
  chain_t *psetj;

  if(accept->m==current->m){
    SAME = TRUE;
    for(k=0; k<accept->m; k++) if(current->tau[k] != accept->tau[k]) SAME = FALSE;
    if(SAME == TRUE){
      alpha = 1;
      profile(FALSE, fnname);
      return alpha;
    }
  }


  psetj = EvaluateProposal(PROPCACHE, current, minseglen, N, g1, g2);
  lZ0 = log(psetj->last->prob) + psetj->first->value;

  psetj = EvaluateProposal(PROPCACHE, accept, minseglen, N, g1, g2);
  lZ1 = log(psetj->last->prob) + psetj->first->value;

  if(current->m == 1){
    if(accept->j == 1){
      accept->j = 0;
    }else{
      accept->j = 1;
    }
    psetj = EvaluateProposal(PROPCACHE, accept, minseglen, N, g1, g2);
    lZ2 = log(psetj->last->prob) + psetj->first->value;  //Z[kjp]
    lmax = (lZ2>lZ1) ? lZ2 : lZ1;
    alpha = lZ0 - lmax - 2*log(2) + log(exp(lmax - lZ1) + exp(lmax - lZ2));
  }else if(accept->m == 1){
    if(current->j == 1){
      current->j = 0;
    }else{
      current->j = 1;
    }
    psetj = EvaluateProposal(PROPCACHE, current, minseglen, N, g1, g2);
    lZ2 = log(psetj->last->prob) + psetj->first->value;  //Z[kj]
    lmax = (lZ2>lZ0) ? lZ2 : lZ0;
    alpha = lmax - lZ1 + 2*log(2) - log(exp(lmax - lZ0) + exp(lmax - lZ2));
  }else{
      alpha = lZ0 - lZ1 + log(current->m) - log(accept->m);
  }
  
  alpha = (alpha>0) ? 1 : exp(alpha);

  profile(FALSE, fnname);
  return alpha;
}


void single_step4(list_t *PROPCACHE, MCMCitem_t *current, 
                 double *g1, struct G2 *g2, int *minseglen, int *N){
  int fnname = 307;
  profile(TRUE, fnname);
  
  chain_t *proposals;
  MCMCitem_t *accept;
  double alpha;
  
  current->j = floor(runif(0,current->m));
  proposals = EvaluateProposal(PROPCACHE, current, minseglen, N, g1, g2);
  accept = sampleProposal(proposals);

  alpha = AcceptanceProb(current, accept, PROPCACHE, g1, g2, minseglen, N);

  if(runif(0,1) < alpha){
    my_free(current->tau);
    current->m = accept->m;
    current->value = accept->value;
    current->j = accept->j;
    current->tau = (int*)my_calloc(current->m,sizeof(int));
    int i;
    for(i=0;i<current->m;i++) current->tau[i] = accept->tau[i];
  }
  DeleteMCMCitem(accept);

  profile(FALSE, fnname);
  return;
}

void Progress(int tick, int max){
  
  int fnname = 309;
  profile(TRUE, fnname);
  
  if(tick==0){
    Rprintf("|");
    if(max == 1) Rprintf("|\n");
    profile(FALSE, fnname);
    return;
  }
  
  double size = ((double) max)/10;
  double this = (double) tick;
  if(floor(this/size) != floor((this+1)/size)) Rprintf("-");
  if(tick==(max-1)) Rprintf("|\n");
  profile(FALSE, fnname);
  return;
}


void batch_iteration(chain_t *CHAIN, int *steps, double *g1, struct G2 *g2, int *minseglen, int *N, 
                     int *progress, int *CACHEMAX){
  int fnname = 310;
  profile(TRUE, fnname);
  
  int i;
  MCMCitem_t *current;
  list_t *PROPCACHE;
  PROPCACHE = makeList(*CACHEMAX);
  for(i=0;i<*steps;i++){
    if(*progress == TRUE) Progress(i, *steps);
    current = CopyMCMCitem(CHAIN->last); //copy last sample
    single_step4(PROPCACHE, current, g1, g2, minseglen, N);
    PushToChain(CHAIN, current);
  }

  DeleteList(PROPCACHE);
  profile(FALSE, fnname);
  return;
}

//--Convergence test

void Tabulate(double *tabm, double *tabt, chain_t *chain, int *N, double *n){
  int fnname = 311;
  profile(TRUE, fnname);
  
  MCMCitem_t *draw;
  draw = chain->first;
  int j;
  n[0] = 0.0;
  n[1] = 0.0;
  while(draw != NULL){
    tabm[draw->m - 1]++;
    n[0]++;
    if(draw->m == 1){
      tabt[*N]++;
      n[1]++;
    }else{
      for(j=0;j<draw->m;j++){
        tabt[draw->tau[j] - 1]++;
        n[1]++;
      }
    }
    draw = draw->next;
  }
  
  profile(FALSE, fnname);
  return;
}

int Compare(double *tab1, double *tab2, double n1, double n2, int len, double *typeI){
  int fnname = 312;
  profile(TRUE, fnname);
  
  double DM, CM;
  int i;

//  int j;
//  for(j=0;j<len;j++){
//    Rprintf("%.0f\t",tab1[j]);
//  }
//  Rprintf("\n");
//  for(j=0;j<len;j++){
//    Rprintf("%.0f\t",tab2[j]);
//  }
//  Rprintf("\n");


  for(i=1;i<len;i++){
    tab1[i] += tab1[i-1];
    tab2[i] += tab2[i-1];
  }
  DM=0.0;
  double tmp;
  for(i=0; i<len; i++){
    tmp = ((tab1[i]/tab1[len-1]) - (tab2[i]/tab2[len-1]));
    if(tmp < 0) tmp *= -1;
    if(tmp > DM) DM = tmp;
  }
  CM = sqrt(-0.5 * log(*typeI/2) * (n1 + n2)/(n1 * n2));
//  Rprintf("DM:%.6f  CM:%.6f\n",DM,CM);

  if(DM < CM){
    profile(FALSE, fnname);
    return TRUE;
  }else{
    profile(FALSE, fnname);
    return FALSE;
  }
}

int ConvTest(chain_t *chain1, chain_t *chain2, double *typeI, int *N, int *maxM){
  int fnname = 313;
  profile(TRUE, fnname);
  
  double *tabm1, *tabm2, *tabt1, *tabt2, *n1, *n2;
  tabm1 = (double*)my_calloc(*maxM,sizeof(double));
  tabm2 = (double*)my_calloc(*maxM,sizeof(double));
  tabt1 = (double*)my_calloc(*N+1,sizeof(double));
  tabt2 = (double*)my_calloc(*N+1,sizeof(double));
  n1 = (double*)my_calloc(2,sizeof(double));
  n2 = (double*)my_calloc(2,sizeof(double));
  Tabulate(tabm1, tabt1, chain1, N, n1);
  Tabulate(tabm2, tabt2, chain2, N, n2);
//  Rprintf("--Mtest--\n");
  int Mtest = Compare(tabm1, tabm2, n1[0], n2[0], *maxM, typeI);
//  Rprintf("--Ttest--\n");
  int Ttest = Compare(tabt1, tabt2, n1[1], n2[1], *N + 1, typeI);
  my_free(tabm1);
  my_free(tabm2);
  my_free(tabt1);
  my_free(tabt2);
  my_free(n1);
  my_free(n2); 
  if((Mtest==TRUE) & (Ttest==TRUE)){
    profile(FALSE, fnname);
    return TRUE;
  }else{
    profile(FALSE, fnname);
    return FALSE;
  }
}

