#ifndef FILE_MODEPCPT
#define FILE_MODEPCPT

#include "PeriodCPT_General.h"
#include "PeriodCPT_MCMCgeneric.h"
#include "PeriodCPT_MCMC.h"

void List_unique_samples(int *chains, int *maxdrawlen, int *ndraws,
                         int *blank, chain_t *unique){
  MCMCitem_t *search = NULL;
  MCMCitem_t *current = NULL;
  for(int i = 0; i < *ndraws; i++){
    int m = 0;
    for(int j=0; j<*maxdrawlen; j++){
      if(chains[i * *maxdrawlen + j] != *blank) m++;
    }
    current = Make_blank_MCMCitem(m, &(chains[i * *maxdrawlen]));
    current->j    = i;   //position in chains
    current->prob = 1.0; //tally data
    search = Find_in_Chain(unique, current);
    if(search == NULL){
      //does not exist in list, so add a copy to chain
      Push_to_chain(unique, Copy_MCMCitem(current));
    }else{
      //iterate counter for case in table
      search->prob += 1.0;
    }
    Delete_MCMCitem(current);
  }
  return;
}




void Mode_pcpt(int *chains, int *maxM, int *ndraws,
              int *blank, int *N, double *g1, segval_t *g2,
              int *mode){

  //Create table of unique pcpt samples, using the prob slot at counter
  chain_t *unique = Make_Chain();
  MCMCitem_t *search = NULL;
  MCMCitem_t *current = NULL;
  List_unique_samples(chains, maxM, ndraws, blank, unique);
  current = unique->first;
  while(current != NULL){
    current->value = Eval_at_CPT(current->tau, current->m, N, g1, g2);
    current = current->next;
  }


  //Searching for mode pcpt
  search = unique->first;
  double maxcount = search->prob;
  current = unique->first->next;
  while(current != NULL){
    if(current->prob > maxcount){
      //Find mode based on count
      maxcount = current->prob;
      search = current;
    }else if(current->prob == maxcount){
      //If two samples have the same count,
      if(current->value > search->value){
        //  first, select the one that has the highest value
        search = current;
      }else if(current->m < search->m){
        //  second, select the one with fewer pcpts
        search = current;
      }else{
        //  third, select based on numeric order of pcpts
        for(int j = 0; (j < current->m) & (search != current); j++){
          if(current->tau[j] > search->tau[j]){
            search = current;
          }
        }
      }
    }
    current = current->next;
  }


  //Pass out mode pcpt (at search) and pad with blanks
  for(int j = 0; j < *maxM; j++) mode[j] = *blank;
  for(int j = 0; j < search->m; j++) mode[j] = search->tau[j];

  //clean-up
  Delete_Chain(unique);
  return;
}


void Calc_Seg_SumStats(MCMCitem_t *mcmc, int seg, double **Stats,
                       int *N, double *SegStats){
  int nStats = (int)Stats[0][0];
  for(int i = 0; i < nStats; i++) SegStats[i] = 0.0;
  int starttau, thistau;
  thistau = mcmc->tau[seg];
  if(seg == 0){
    starttau = mcmc->tau[mcmc->m - 1];
  }else{
    starttau = mcmc->tau[seg - 1];
  }
  //starttau++;
  //if(starttau > *N) starttau -= *N;

  if(mcmc->m == 1){ //Do first outside loop to avoid failing while condition for first pcpt!
    for(int k = 0; k<nStats; k++){
      SegStats[k] += Stats[k+1][thistau-1];
    }
    thistau--;
    if(thistau <= 0) thistau += *N;
  }

  while(thistau != starttau){
    for(int k = 0; k<nStats; k++){
      SegStats[k] += Stats[k+1][thistau-1];
    }
    thistau--;
    if(thistau <= 0) thistau += *N;
  }


  return;
}


void Evaluate_fits(MCMCitem_t *draw, char **Pdist, double *Phyp,
                   double **sumStats, int *N, int *nfits, double *fits, int *err){
  //Get Fit_FN for specified sampling distribution
  Fit_FN_Ptr *Fit_FN;
  Get_Fit_FN(Pdist, &Fit_FN, err);
  if(*err != 0) return;

  //Initialise fits output
  for(int i = 0; i < *nfits; i++) fits[i] = 0.0;


  int nStats = (int)sumStats[0][0];
  for(int j = 0; j < draw->m; j++){
    //Evaluate summary stats for segment
    double *thisStats = my_calloc(nStats, sizeof(double));
    Calc_Seg_SumStats(draw, j, sumStats, N, thisStats);
    //Calc fit on segment, add into fits output
    Fit_FN(thisStats, nStats, N, Phyp, err, fits);
    my_free(thisStats);
  }

  for(int i = 0; i < *nfits; i++) fits[i] *= -2;
  return;
}


void Evaluate_SegParamMode(MCMCitem_t *mcmc, char **Pdist, double *Phyp, int *N,
                           int *maxM, int *blank, double **sumStats, int *nsegparams,
                           double *params, int *err){

  Param_Mode_Ptr *Param_Mode;
  Get_Param_Mode(Pdist, &Param_Mode, err);
  if(*err != 0) return;

  for(int i=0; i<(*maxM * *nsegparams); i++) params[i] = *blank;

  for(int seg = 0; seg < mcmc->m; seg++){
    double *thisStats = my_calloc(sumStats[0][0],sizeof(double));
    Calc_Seg_SumStats(mcmc, seg, sumStats, N, thisStats);
    Param_Mode(thisStats, sumStats[0][0], Phyp, &(params[seg * *nsegparams]), err);
    my_free(thisStats);
  }
  return;
}



void Mode_Fit_Calcuation(double *data, int *time, int *n, int *N, char **Pdist, double *Phyp,
                      int *tau, int *blank, int *maxM, int *err, double *fits, int *nfits,
                      double *params, int *nsegparams){

  Summary_Stats_Ptr *Summary_Stats;
  Get_SumStats_FN(Pdist, &Summary_Stats, err);
  if(*err != 0) return;
  double **SumStats = Summary_Stats(data, time, n, N);

  MCMCitem_t *draw = my_calloc(1,sizeof(MCMCitem_t));
  draw->m = 0;
  for(int i=0; i<*maxM;i++){
    if(tau[i] != *blank) draw->m++;
  }
  draw->tau = my_calloc(draw->m, sizeof(int));
  for(int i=0; i<draw->m; i++) draw->tau[i] = tau[i];

  Evaluate_SegParamMode(draw, Pdist, Phyp, N, maxM, blank, SumStats, nsegparams, params, err);
  Evaluate_fits(draw, Pdist, Phyp, SumStats, N, nfits, fits, err);

  Delete_MCMCitem(draw);

}




#endif //FILE_MODEPCPT
