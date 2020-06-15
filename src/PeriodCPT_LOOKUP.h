#ifndef FILE_LOOKUP
#define FILE_LOOKUP

#include "PeriodCPT_General.h"
#include "PeriodCPT_distributions.h"
#include "PeriodCPT_CaseFNs.h"

void Evaluate_G1(Mprior_Ptr Mprior, int *maxM, int *N, double *spread,
                 double *Mhyp, int *minseglen, double *g1){
  //Rprintf("Evaluate_G1\n");
  for(int m = 1; m <= *maxM; m++){
    //model prior
    g1[m-1] = Mprior(m, maxM, Mhyp);	//Prior for number of segments
    if(m > 1){
      g1[m-1] -= log((double) *N);  		//Prior for anchor cpt(?)
      //Global contribution to the excess spacing between segments
      g1[m-1] += lgammafn( (double) (*N - (m * *minseglen) + 1) );
      g1[m-1] += lgammafn( (double)m * *spread);
      g1[m-1] -= lgammafn( (double)(*N - (m * *minseglen)) +
                 ((double)m * *spread));
    }
  }
  return;
}

void Evaluate_G2(double **SumStats, Samp_Dist_Ptr Samp_Dist, int *N,
                 int *minseglen, double *spread, double *Phyp,
                 segval_t *g2){
  //Rprintf("Evaluate_G2\n");
  //Evaluate lookup table for function g2(tau, len)
  int id, startSeg;
  double *Spart     = my_calloc(*N, sizeof(double));
  int nStats        = (int)SumStats[0][0];
  double *thisStats = my_calloc(nStats, sizeof(double));

  id = 0;
  for(int tau = 1; tau <= *N; tau++){
    for(int s = 0; s < nStats; s++){
      thisStats[s] = 0.0;
    }

    for(int len = 1; len <= *N; len++){

      if(tau == 1){ //calc & store the per segment contribution to
                    //  the excess spacing between segments
        if((len < *N) & (len >= *minseglen)){
          Spart[len-1] = lgammafn((double)(len - *minseglen) + *spread) ;
          Spart[len-1] -= lgammafn((double)(len - *minseglen) + 1);
          Spart[len-1] -= lgammafn( *spread );
        }else{
          Spart[len-1] = 0.0;
        }
      }

      //Evaluate summary stats for segment (i.e. update last
      //   calc with extra obs at startSeg)
      startSeg = tau - len + 1;
      if(startSeg<=0) startSeg += *N;   //range: 1--N
      for(int s = 0; s < nStats; s++){
        thisStats[s] += SumStats[s+1][startSeg - 1];
      }

      if(len < *minseglen){  //Immediately deal with lengths shorter
                             //  than minseglen
        g2[id].cpt    = tau;
        g2[id].seglen = len;
        g2[id].value  = -INFINITY;
      }else{
        g2[id].cpt    = tau;
        g2[id].seglen = len;
        g2[id].value  = Spart[len-1] + Samp_Dist(thisStats, nStats, Phyp);
      }
      id++;
    }
  }
  my_free(Spart);
  my_free(thisStats);
  return;
}

void MAKE_LOOK_TABLES(
    double    *data,         //Periodic time series data
    int       *time,         //Within period time index
    int       *n,            //Length of data
    int       *N,            //Period length
    int       *minseglen,    //Minimum segment length
    int       *maxM,         //max number of possible pcpts: floor(N/l)
    char     **Mdist,        //character of pcpt total prior distribution
    double    *Mhyp,         //Hyper-parameter(s) for pcpt total prior
    double    *spread,       //Hyper-parameter for pcpt location prior
    char     **Pdist,        //function of sampling distribution
    double    *Phyp,         //Hyper-para for segment parameter priors
    int       *err,          //Error flag
    double    *g1,           //Lookup list (m)
    segval_t  *g2){          //Lookup table (pcpt)
  //Rprintf("MAKE_LOOK_TABLES\n");
  Mprior_Ptr *Mprior = NULL;
  Samp_Dist_Ptr *Samp_Dist = NULL;
  Summary_Stats_Ptr *Summary_Stats = NULL;
  Get_Mprior(Mdist, &Mprior, err);
  Get_Pprior(Pdist, &Samp_Dist, err);
  Get_SumStats_FN(Pdist, &Summary_Stats, err);

  if(*err != 0){
    return;
  }
  Evaluate_G1(Mprior, maxM, N, spread, Mhyp, minseglen, g1);
  double **SumStats = Summary_Stats(data, time, n, N);
  Evaluate_G2(SumStats, Samp_Dist, N, minseglen, spread, Phyp, g2);
  Free_Summary_Stats(SumStats);

  return;
}


double Get_g2(int tau, int len, segval_t *g2, int *N){
  //Rprintf("Get_g2\n");
  //order: t1l1 t1l2 t1l3 t2l1 t2l2 t2l3 t3l1 t3l2 t3l3
  //Includes -Inf cases where len<minseglen (retain for spacing!!!)
  return g2[(tau-1) * *N + len - 1].value;
}

double Get_g1(int m, double *g1){
  //Rprintf("Get_g1\n");
  return g1[m-1];
}

double Eval_at_CPT(int *tau, int m, int *N, double *g1, segval_t *g2){
  //Rprintf("Eval_at_CPT\n");
  double out = Get_g1(m,g1);
  if(m > 1){
    out += Get_g2(tau[0], tau[0] - tau[m-1] + *N, g2, N);
    for(int j = 1; j < m; j++){
      out += Get_g2(tau[j], tau[j] - tau[j-1], g2, N);
    }
  }else{
    out += Get_g2(tau[0], *N, g2, N);
  }
  return out;
}

/*
void PrintLookup(double *g1, segval_t *g2, int *maxM, int *N){
  //Rprintf("PrintLookup\n");
  FILE *f;

  //---------------------------------------

  f = fopen("g1table.txt","w");
  fprintf(f, "model, g1\n");
  for(int id = 0; id < *maxM; id++){
    fprintf(f, "%d, %f\n", id+1, g1[id]);
  }
  fclose(f);

  //---------------------------------------

  f = fopen("g2table.txt","w");
  fprintf(f,"tau, seglen, value\n");
  for(int id = 0; id < (*N * *N); id++){
    fprintf(f, "%d, %d, %f\n", g2[id].cpt, g2[id].seglen, g2[id].value);
  }
  fclose(f);

  //---------------------------------------
  return;
}
*/


#endif //FILE_LOOKUP
