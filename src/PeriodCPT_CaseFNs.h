#ifndef FILE_CASEFNS
#define FILE_CASEFNS

#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <math.h>

#include "PeriodCPT_distributions.h"
#include "PeriodCPT_General.h"

//-------------------------------------------------------------
// Calloc and free functions for summary statistics function
//--------------------------------------------------------------

void Free_Summary_Stats(double **SumStats){
  //Rprintf("Free_Summary_Stats\n");
  for(int i = 1; i <= SumStats[0][0]; i++){
    my_free(SumStats[i]);
  }
  my_free(SumStats[0]);
  my_free(SumStats);
  return;
}

double **Make_Summary_Stats(int number_of_stat_types, int length){
  //Rprintf("Make_Summary_Stats\n");
  double **SumStats, *tmp;
  SumStats = (double **)my_calloc(1 + number_of_stat_types, sizeof(double *));
  tmp = (double *)my_calloc(1, sizeof(double));
  tmp[0] = (double)number_of_stat_types;
  SumStats[0] = tmp;
  for(int i = 1; i <= number_of_stat_types; i++){
    SumStats[i] = (double *)my_calloc(2 * length, sizeof(double));
  }
  return SumStats;
}

//------------------------------------------------------------
// Prior functions for the number of periodic changepoint
//------------------------------------------------------------

double Mprior_pois(int m, int *max, double *params){
  //Rprintf("Mprior_pois\n");
  //Trunc-Pois on interval 1, ..., max.
  double rate = params[0];
  return dTpois(m, *max, rate, TRUE);
}

double Mprior_unif(int m, int *max, double *params){
  //Rprintf("Mprior_unif\n");
  //Discrete-Unif on interval 1, ..., max.
  return -log(*max);
}

double Mprior_BLANK(int m, int *max, double *params){
  //Rprintf("Mprior_BLANK\n");
  //Mprior blank holding funciton for errors.
  return 0.0;
}

//------------------------------------------------------------
//Functions related to the sampling distribution
//------------------------------------------------------------

//-BLANK------------------------------------------------------

double Samp_Dist_BLANK(double *SumStats, int nStat, double *Phyp){
  //Rprintf("Samp_Dist_BLANK\n");
  //Samp_Dist blank holding funciton for errors.
  return 0.0;
}

double **Summary_Stats_BLANK(double *data, int *time, int *n, int *N){
  //Rprintf("Summary_Stats_BLANK\n");
  //Summary_Stats blank holding funciton for errors.
  return NULL;
};

void Sufficient_Stats_BLANK(double *SumStats, int nSumm, int nSuff, double *Phyp, double *out){
  //Rprintf("Sufficient_Stats_BLANK\n");
  return;
}

void Param_Mode_BLANK(double *SumStats, int nSumm, double *Phyp, double* Pmode, int *err){
  return;
}

void Fit_FN_BLANK(double *SumStats, int nSumm, int *N, double *Phyp, int *err, double *fits){
  double *Pmode = (double *)calloc(1, sizeof(double));
  Param_Mode_BLANK(SumStats, nSumm, Phyp, Pmode, err);

  //Function used in MCMC:
  //log\{ \int f(y|\theta,\tau) \pi(\theta|\tau) d\theta \}
  fits[0] += Samp_Dist_BLANK(SumStats, nSumm, Phyp);

  //maximised (log-)joint of segment likelihood and segment prior (ie evaluated at posterior mode):
  //\max_\theta [ log\{f(y|\theta,\tau) \pi(\theta|\tau) \} ]
  fits[1] += 0.0;

  //(log-)integrated segment sampling distribution:
  // \log\{ \int f(y|\theta, \tau) d\theta \}
  fits[2] += 0.0;

  //(log-)segment likelihood evaluated at posterior mode:
  // \log\{ f(y|\theta^*, \tau) \} where \theta^* = \argmax{\pi(\theta|y,\tau)}
  fits[3] += 0.0;

  //maximised (log-)segment likelihood evaluated:
  fits[4] += 0.0;

  free(Pmode);
  return;
}













//-bern------------------------------------------------------

double Samp_Dist_bern(double *SumStats, int nStat, double *Phyp){
  //Rprintf("Samp_Dist_bern\n");
  double sum_x          = SumStats[0];
  double count_x        = SumStats[1];
  double alpha          = Phyp[0];
  double beta           = Phyp[1];
  return dbbinom(sum_x, count_x, alpha, beta, TRUE) +
    dmhg(1.0, sum_x, count_x, TRUE);
}

double **Summary_Stats_bern(double *data, int *time, int *n, int *N){
  //Rprintf("Summary_Stats_bern\n");
  double **SumStats = Make_Summary_Stats(2, *N);
  for(int i = 0; i < *n; i++){
    SumStats[1][time[i] - 1] += data[i];    //0: sum_xi
    SumStats[2][time[i] - 1] += 1;          //1: count_xi
  }
  return SumStats;
}

void Sufficient_Stats_bern(double *SumStats, int nSumm, int nSuff, double *Phyp, double *out){
  //Rprintf("Sufficient_Stats_bern\n");
  double sum_x          = SumStats[0];
  double count_x        = SumStats[1];
  double alpha          = Phyp[0];
  double beta           = Phyp[1];
  out[0] = alpha + sum_x;
  out[1] = beta  + count_x - sum_x;
  return;
}

void Param_Mode_bern(double *SumStats, int nSumm, double *Phyp, double* Pmode, int *err){
  double *SuffStats = (double *)calloc(2, sizeof(double));
  Sufficient_Stats_bern(SumStats, nSumm, 2, Phyp, SuffStats);
  if( SuffStats[0]>1 & SuffStats[1]>1 ){
    Pmode[0] = (SuffStats[0]-1)/(SuffStats[0] + SuffStats[1]-2);
  }else if(SuffStats[0] <= 1 & SuffStats[1]>1){
    Pmode[0] = 0;
  }else if(SuffStats[0] > 1 & SuffStats[1]<=1){
    Pmode[0] = 1;
  }else{
    Pmode[0] = 0;
    *err = 101; //Mode is not unique! Dist is either flat or bimodal at {0,1}.
                //Consider changing Prior hyper-params to be > 1.
  }
  free(SuffStats);
  return;
}

void Fit_FN_bern(double *SumStats, int nSumm, int *N, double *Phyp, int *err, double *fits){
  double sum_x          = SumStats[0];
  double count_x        = SumStats[1];
  double alpha          = Phyp[0];
  double beta           = Phyp[1];
  double *Pmode = (double *)calloc(1, sizeof(double));
  Param_Mode_bern(SumStats, nSumm, Phyp, Pmode, err);

  //Function used in MCMC:
  //log\{ \int f(y|\theta,\tau) \pi(\theta|\tau) d\theta \}
  fits[0] += Samp_Dist_bern(SumStats, nSumm, Phyp);

  //maximised (log-)joint of segment likelihood and segment prior (ie evaluated at posterior mode):
  //\max_\theta [ log\{f(y|\theta,\tau) \pi(\theta|\tau) \} ]
  fits[1] += (sum_x + Phyp[0] - 1) * log(Pmode[0]) + (count_x - sum_x + Phyp[1] - 1) * log(1-Pmode[0]) +
    lgamma(Phyp[0] + Phyp[1]) - lgammafn(Phyp[0]) - lgammafn(Phyp[1]);

  //(log-)integrated segment sampling distribution:
  // \log\{ \int f(y|\theta, \tau) d\theta \}
  fits[2] += lgammafn(sum_x + 1) + lgammafn(count_x - sum_x + 1) - lgammafn(count_x + 2);

  //(log-)segment likelihood evaluated at posterior mode:
  // \log\{ f(y|\theta^*, \tau) \} where \theta^* = \argmax{\pi(\theta|y,\tau)}
  fits[3] += sum_x * log(Pmode[0]) + (count_x - sum_x) * log(1 - Pmode[0]);

  //maximised (log-)segment likelihood evaluated:
  fits[4] += sum_x * log(sum_x/count_x) + (count_x - sum_x) * log(1 - sum_x/count_x);


  free(Pmode);
  return;
}












//-binom-----------------------------------------------------


double Samp_Dist_binom(double *SumStats, int nStat, double *Phyp){
  //Rprintf("Samp_Dist_binom\n");
  double sum_x          = SumStats[0];
  double sum_n          = SumStats[1];
  double sum_lchoose_nx = SumStats[2];
  double alpha          = Phyp[0];
  double beta           = Phyp[1];
  return dbbinom(sum_x, sum_n, alpha, beta, TRUE) +
    dmhg(sum_lchoose_nx, sum_x, sum_n, TRUE);
}

double **Summary_Stats_binom(double *data, int *time, int *n, int *N){
  //Rprintf("Summary_Stats_binom\n");
  double **SumStats = Make_Summary_Stats(3, *N);
  for(int i = 0; i < *n; i++){
    SumStats[1][time[i] - 1] += data[i];                        //0: sum_xi
    SumStats[2][time[i] - 1] += data[i + *n];                   //1: sum_ni
    SumStats[3][time[i] - 1] += lchoose(data[i + *n], data[i]); //2: sum_log(choose(ni, xi))
  }
  return SumStats;
}

void Sufficient_Stats_binom(double *SumStats, int nSumm, int nSuff, double *Phyp, double *out){
  //Rprintf("Sufficient_Stats_binom\n");
  double sum_x          = SumStats[0];
  double sum_n          = SumStats[1];
  //double sum_lchoose_nx = SumStats[2];
  double alpha          = Phyp[0];
  double beta           = Phyp[1];
  out[0] = alpha + sum_x;
  out[1] = beta  + sum_n - sum_x;
  return;
}

void Param_Mode_binom(double *SumStats, int nSumm, double *Phyp, double* Pmode, int *err){
  double *SuffStats = (double *)calloc(2, sizeof(double));
  Sufficient_Stats_binom(SumStats, nSumm, 2, Phyp, SuffStats);
  if( SuffStats[0]>1 & SuffStats[1]>1 ){
    Pmode[0] = (SuffStats[0]-1)/(SuffStats[0] + SuffStats[1]-2);
  }else if(SuffStats[0] <= 1 & SuffStats[1]>1){
    Pmode[0] = 0;
  }else if(SuffStats[0] > 1 & SuffStats[1]<=1){
    Pmode[0] = 1;
  }else{
    Pmode[0] = 0;
    *err = 101; //Mode is not unique! Dist is either flat or bimodal at {0,1}.
                //Consider changing Prior hyper-params to be > 1.
  }
  free(SuffStats);
  return;
}

void Fit_FN_binom(double *SumStats, int nSumm, int *N, double *Phyp, int *err, double *fits){
  double sum_x          = SumStats[0];
  double sum_n          = SumStats[1];
  double sum_lchoose_nx = SumStats[2];
  double alpha          = Phyp[0];
  double beta           = Phyp[1];
  double *Pmode = (double *)calloc(1, sizeof(double));
  Param_Mode_binom(SumStats, nSumm, Phyp, Pmode, err);

  //Function used in MCMC:
  //log\{ \int f(y|\theta,\tau) \pi(\theta|\tau) d\theta \}
  fits[0] += Samp_Dist_binom(SumStats, nSumm, Phyp);

  //maximised (log-)joint of segment likelihood and segment prior (ie evaluated at posterior mode):
  //\max_\theta [ log\{f(y|\theta,\tau) \pi(\theta|\tau) \} ]
  fits[1] += (sum_x + Phyp[0] - 1) * log(Pmode[0]) + (sum_n - sum_x + Phyp[1] - 1) * log(1-Pmode[0]) +
    lgamma(Phyp[0] + Phyp[1]) - lgammafn(Phyp[0]) - lgammafn(Phyp[1]) + sum_lchoose_nx;

  //(log-)integrated segment sampling distribution:
  // \log\{ \int f(y|\theta, \tau) d\theta \}
  fits[2] += lgammafn(sum_x + 1) + lgammafn(sum_n - sum_x + 1) - lgammafn(sum_n + 2) + sum_lchoose_nx;

  //(log-)segment likelihood evaluated at posterior mode:
  // \log\{ f(y|\theta^*, \tau) \} where \theta^* = \argmax{\pi(\theta|y,\tau)}
  fits[3] += sum_lchoose_nx + sum_x * log(Pmode[0]) + (sum_n - sum_x) * log(1 - Pmode[0]);

  //maximised (log-)segment likelihood evaluated:
  fits[4] += sum_lchoose_nx + sum_x * log(sum_x/sum_n) + (sum_n - sum_x) * log(1 - sum_x/sum_n);

  free(Pmode);
  return;
}










//-pois-----------------------------------------------------

double Samp_Dist_pois(double *SumStats, int nStat, double *Phyp){
  //Rprintf("Samp_Dist_pois\n");
  double sum_x          = SumStats[0];
  double count_x        = SumStats[1];
  double sum_lfact      = SumStats[2];
  double alpha          = Phyp[0];
  double beta           = Phyp[1];
  return dnbinom(sum_x, alpha, beta/(beta + count_x), TRUE) +
    dmult(count_x, sum_x, sum_lfact, TRUE);
}

double **Summary_Stats_pois(double *data, int *time, int *n, int *N){
  //Rprintf("Summary_Stats_pois\n");
  double **SumStats = Make_Summary_Stats(3, *N);
  for(int i = 0; i < *n; i++){
    SumStats[1][time[i] - 1] += data[i];            //0: sum_xi
    SumStats[2][time[i] - 1] += 1;                  //1: count_xi
    SumStats[3][time[i] - 1] += lgamma(data[i]+1);  //2: sum_log(xi!)
  }
  return SumStats;
}

void Sufficient_Stats_pois(double *SumStats, int nSumm, int nSuff, double *Phyp, double *out){
  //Rprintf("Sufficient_Stats_pois\n");
  double sum_x          = SumStats[0];
  double count_x        = SumStats[1];
  //double sum_lfact      = SumStats[2];
  double alpha          = Phyp[0];
  double beta           = Phyp[1];
  out[0] = alpha + sum_x;
  out[1] = beta  + count_x;
  return;
}

void Param_Mode_pois(double *SumStats, int nSumm, double *Phyp, double* Pmode, int *err){
  double *SuffStats = (double *)calloc(2, sizeof(double));
  Sufficient_Stats_pois(SumStats, nSumm, 2, Phyp, SuffStats);
  if( SuffStats[0]>=1){
    Pmode[0] = (SuffStats[0] - 1)/SuffStats[1];
  }else{
    Pmode[0] = 0;
  }
  free(SuffStats);
  return;
}

void Fit_FN_pois(double *SumStats, int nSumm, int *N, double *Phyp, int *err, double *fits){
  double sum_x          = SumStats[0];
  double count_x        = SumStats[1];
  double sum_lfact      = SumStats[2];
  double alpha          = Phyp[0];
  double beta           = Phyp[1];
  double *Pmode = (double *)calloc(1, sizeof(double));
  Param_Mode_pois(SumStats, nSumm, Phyp, Pmode, err);

  //Function used in MCMC:
  //log\{ \int f(y|\theta,\tau) \pi(\theta|\tau) d\theta \}
  fits[0] += Samp_Dist_pois(SumStats, nSumm, Phyp);

  //maximised (log-)joint of segment likelihood and segment prior (ie evaluated at posterior mode):
  //\max_\theta [ log\{f(y|\theta,\tau) \pi(\theta|\tau) \} ]
  fits[1] += -sum_lfact + Phyp[0]*log(Phyp[1]) - lgammafn(Phyp[0]) + (sum_x + Phyp[0] - 1)*log(Pmode[0]) -
    (Phyp[1] + count_x)*Pmode[0];

  //(log-)integrated segment sampling distribution:
  // \log\{ \int f(y|\theta, \tau) d\theta \}
  fits[2] += -sum_lfact + lgammafn(sum_x + 1) - (sum_x + 1)*log(count_x);

  //(log-)segment likelihood evaluated at posterior mode:
  // \log\{ f(y|\theta^*, \tau) \} where \theta^* = \argmax{\pi(\theta|y,\tau)}
  fits[3] += -sum_lfact + sum_x * log(Pmode[0]) - count_x*Pmode[0];

  //maximised (log-)segment likelihood evaluated:
  fits[4] += -sum_lfact + sum_x * log(sum_x/count_x) - sum_x;

  free(Pmode);
  return;
}








//-mean-----------------------------------------------------

double Samp_Dist_mean(double *SumStats, int nStat, double *Phyp){
  //Rprintf("Samp_Dist_mean\n");
  double count_x        = SumStats[0];
  double sum_x          = SumStats[1];
  double sum_xx         = SumStats[2];
  double m              = Phyp[0];
  double c              = Phyp[1];
  double var            = Phyp[2];
  return dmvnorm(count_x, sum_x, sum_xx, m, 1/c, 1/var, TRUE);
}

double **Summary_Stats_mean(double *data, int *time, int *n, int *N){
  //Rprintf("Summary_Stats_mean\n");
  double **SumStats = Make_Summary_Stats(3, *N);
  for(int i = 0; i < *n; i++){
    SumStats[1][time[i] - 1] += 1;                  //0: count_xi
    SumStats[2][time[i] - 1] += data[i];            //1: sum_xi
    SumStats[3][time[i] - 1] += data[i]*data[i];    //2: sum_xi^2
  }
  return SumStats;
}

void Sufficient_Stats_mean(double *SumStats, int nSumm, int nSuff, double *Phyp, double *out){
  //Rprintf("Sufficient_Stats_mean\n");
  double count_x        = SumStats[0];
  double sum_x          = SumStats[1];
  //double sum_xx         = SumStats[2];
  double m              = Phyp[0];
  double c              = Phyp[1];
  double var            = Phyp[2];
  double tmp = count_x*c + var;
  out[0] = (c*sum_x + var*m)/tmp;
  out[1] = var*c/tmp;
  return;
}

void Param_Mode_mean(double *SumStats, int nSumm, double *Phyp, double* Pmode, int *err){
  double *SuffStats = (double *)calloc(2, sizeof(double));
  Sufficient_Stats_mean(SumStats, nSumm, 2, Phyp, SuffStats);
  Pmode[0] = SumStats[0];
  free(SuffStats);
  return;
}

void Fit_FN_mean(double *SumStats, int nSumm, int *N, double *Phyp, int *err, double *fits){
  double count_x        = SumStats[0];
  double sum_x          = SumStats[1];
  double sum_xx         = SumStats[2];
  double m              = Phyp[0];
  double c              = Phyp[1];
  double var            = Phyp[2];
  double *Pmode = (double *)calloc(1, sizeof(double));
  Param_Mode_mean(SumStats, nSumm, Phyp, Pmode, err);

  //Function used in MCMC:
  //log\{ \int f(y|\theta,\tau) \pi(\theta|\tau) d\theta \}
  fits[0] += Samp_Dist_mean(SumStats, nSumm, Phyp);

  //maximised (log-)joint of segment likelihood and segment prior (ie evaluated at posterior mode):
  //\max_\theta [ log\{f(y|\theta,\tau) \pi(\theta|\tau) \} ]
  fits[1] += 0.5*count_x*log(2*PI*var) - 0.5*(sum_xx - 2*sum_x*Pmode[0] + count_x*Pmode[0]*Pmode[0])/var +
    0.5*log(2*PI*c) - (Pmode[0] - m)*(Pmode[0] - m)*0.5/c;

  //(log-)integrated segment sampling distribution:
  // \log\{ \int f(y|\theta, \tau) d\theta \}
  fits[2] += -0.5*log(count_x) - 0.5*(count_x-1)*log(2*PI*var) - 0.5*(sum_xx - sum_x*(sum_x/count_x))/var;

  //(log-)segment likelihood evaluated at posterior mode:
  // \log\{ f(y|\theta^*, \tau) \} where \theta^* = \argmax{\pi(\theta|y,\tau)}
  fits[3] += 0.5*count_x*log(2*PI*var) - 0.5*(sum_xx - 2*sum_x*Pmode[0] + count_x*Pmode[0]*Pmode[0])/var;

  //maximised (log-)segment likelihood evaluated:
  fits[4] += 0.5*count_x*log(2*PI*var) - 0.5*(sum_xx - sum_x*(sum_x/count_x))/var;

  free(Pmode);
  return;
}








//-norm-----------------------------------------------------

double Samp_Dist_norm(double *SumStats, int nStat, double *Phyp){
  //Rprintf("Samp_Dist_norm\n");
  double count_x = SumStats[0];
  double sum_x   = SumStats[1];
  double sum_xx  = SumStats[2];
  double m       = Phyp[0];
  double c       = Phyp[1];
  double alpha   = Phyp[2];
  double beta    = Phyp[3];
  return dmvstMeanVar(count_x, sum_x, sum_xx, m, 1/c, alpha, beta, TRUE);
}

double ** Summary_Stats_norm(double *data, int *time, int *n, int *N){
  //Rprintf("Summary_Stats_norm\n");
  double **SumStats = Make_Summary_Stats(3, *N);
  for(int i = 0; i < *n; i++){
    SumStats[1][time[i] - 1] += 1;                  //0: count_xi
    SumStats[2][time[i] - 1] += data[i];            //1: sum_xi
    SumStats[3][time[i] - 1] += data[i]*data[i];    //2: sum_xi^2
  }
  return SumStats;
}

void Sufficient_Stats_norm(double *SumStats, int nSumm, int nSuff, double *Phyp, double *out){
  //Rprintf("Sufficient_Stats_norm\n");
  double count_x = SumStats[0];
  double sum_x   = SumStats[1];
  double sum_xx  = SumStats[2];
  double m       = Phyp[0];
  double c       = Phyp[1];
  double alpha   = Phyp[2];
  double beta    = Phyp[3];
  double tmp = count_x * c + 1;
  out[0] = (c*sum_x + m) / tmp;
  out[1] = c/tmp;
  out[2] = alpha + 0.5 * count_x;
  out[3] = beta + 0.5 * (sum_xx + sum_x * (sum_x/count_x)) +
    0.5 * count_x * (sum_x/count_x + m) * (sum_x/count_x + m) / tmp;
  return;
}

void Param_Mode_norm(double *SumStats, int nSumm, double *Phyp, double* Pmode, int *err){
  double *SuffStats = (double *)calloc(2, sizeof(double));
  Sufficient_Stats_mean(SumStats, nSumm, 2, Phyp, SuffStats);
  Pmode[0] = SumStats[0];
  Pmode[1] = SumStats[3]/(SumStats[2] + 1.5);
  free(SuffStats);
  return;
}

void Fit_FN_norm(double *SumStats, int nSumm, int *N, double *Phyp, int *err, double *fits){
  double count_x = SumStats[0];
  double sum_x   = SumStats[1];
  double sum_xx  = SumStats[2];
  double m       = Phyp[0];
  double c       = Phyp[1];
  double alpha   = Phyp[2];
  double beta    = Phyp[3];
  double *Pmode = (double *)calloc(2, sizeof(double));
  Param_Mode_norm(SumStats, nSumm, Phyp, Pmode, err);

  //Function used in MCMC:
  //log\{ \int f(y|\theta,\tau) \pi(\theta|\tau) d\theta \}
  fits[0] += Samp_Dist_norm(SumStats, nSumm, Phyp);

  //maximised (log-)joint of segment likelihood and segment prior (ie evaluated at posterior mode):
  //\max_\theta [ log\{f(y|\theta,\tau) \pi(\theta|\tau) \} ]
  fits[1] += alpha*log(beta) - lgammafn(alpha) - (alpha + 1) * log(Pmode[1]) - beta/Pmode[1] +
    -0.5*log(2*PI*c*Pmode[1]) - 0.5*(Pmode[0]-m)*(Pmode[0]-m) +
    -0.5*count_x*log(2*PI*Pmode[1]) - 0.5*(sum_xx - 2*Pmode[0]*sum_x + count_x*Pmode[0]*Pmode[0])/Pmode[1];

  //(log-)integrated segment sampling distribution:
  // \log\{ \int f(y|\theta, \tau) d\theta \}
  if(count_x < 3){
    fits[2] += log(0);
    *err = 102;  //Small sample size for reliable segment estimate.
                 //  Consider increasing minimum segment length

  }else{
    fits[2] += -0.5*(count_x-1)*log(2*PI) - 0.5*log(count_x) + lgammafn(0.5*(count_x-3)) -
      (0.5*(count_x-3)) * log(0.5*(sum_xx - sum_x*(sum_x/count_x)));
  }

  //(log-)segment likelihood evaluated at posterior mode:
  // \log\{ f(y|\theta^*, \tau) \} where \theta^* = \argmax{\pi(\theta|y,\tau)}
  fits[3] += -0.5*count_x*log(2*PI*Pmode[1]) - 0.5*(sum_xx - 2*Pmode[0]*sum_x +
    count_x*Pmode[0]*Pmode[0])/Pmode[1];

  //maximised (log-)segment likelihood evaluated:
  fits[5] += -0.5*count_x - 0.5*count_x*log(2 * PI *(sum_xx/count_x - (sum_x/count_x)*(sum_x/count_x)));

  free(Pmode);
  return;
}















//-var-----------------------------------------------------

double Samp_Dist_var(double *SumStats, int nStat, double *Phyp){
  //Rprintf("Samp_Dist_var\n");
  double count_x = SumStats[0];
  double sum_xx  = SumStats[1];
  double alpha   = Phyp[2];
  double beta    = Phyp[3];
  return dmvstvar(count_x, 0.0, sum_xx, 0.0, alpha, beta, TRUE);
}

double ** Summary_Stats_var(double *data, int *time, int *n, int *N){
  //Rprintf("Summary_Stats_var\n");
  double **SumStats = Make_Summary_Stats(3, *N);
  for(int i = 0; i < *n; i++){
    SumStats[1][time[i] - 1] += 1;                  //0: count_xi
    SumStats[2][time[i] - 1] += data[i]*data[i];    //1: sum_xi^2
  }
  return SumStats;
}

void Sufficient_Stats_var(double *SumStats, int nSumm, int nSuff, double *Phyp, double *out){
  //Rprintf("Sufficient_Stats_var\n");
  double count_x = SumStats[0];
  double sum_xx  = SumStats[1];
  double alpha   = Phyp[2];
  double beta    = Phyp[3];
  out[0] = alpha + 0.5 * count_x;
  out[1] = beta + 0.5 * sum_xx;
  return;
}

void Param_Mode_var(double *SumStats, int nSumm, double *Phyp, double* Pmode, int *err){
  double *SuffStats = (double *)calloc(2, sizeof(double));
  Sufficient_Stats_mean(SumStats, nSumm, 2, Phyp, SuffStats);
  Pmode[0] = SuffStats[1]/(SuffStats[0]+1);
  free(SuffStats);
  return;
}

void Fit_FN_var(double *SumStats, int nSumm, int *N, double *Phyp, int *err, double *fits){
  double count_x = SumStats[0];
  double sum_xx  = SumStats[1];
  double alpha   = Phyp[2];
  double beta    = Phyp[3];
  double *Pmode = (double *)calloc(1, sizeof(double));
  Param_Mode_var(SumStats, nSumm, Phyp, Pmode, err);

  //Function used in MCMC:
  //log\{ \int f(y|\theta,\tau) \pi(\theta|\tau) d\theta \}
  fits[0] += Samp_Dist_var(SumStats, nSumm, Phyp);

  //maximised (log-)joint of segment likelihood and segment prior (ie evaluated at posterior mode):
  //\max_\theta [ log\{f(y|\theta,\tau) \pi(\theta|\tau) \} ]
  fits[1] += alpha*log(beta) - lgammaf(alpha) - (alpha + 1)*log(Pmode[0]) - beta/Pmode[0] +
    -0.5*count_x*log(2*PI*Pmode[0]) - 0.5*sum_xx/Pmode[0];

  //(log-)integrated segment sampling distribution:
  // \log\{ \int f(y|\theta, \tau) d\theta \}
  if(count_x < 2){
    fits[2] += log(0);
    *err = 102;  //Small sample size for reliable segment estimate.
                 //  Consider increasing minimum segment length
  }else{
    fits[2] += -0.5*count_x*log(2*PI) + lgammafn(0.5*(count_x-2)) - 0.5*(count_x-2)*log(0.5*sum_xx);
  }

  //(log-)segment likelihood evaluated at posterior mode:
  // \log\{ f(y|\theta^*, \tau) \} where \theta^* = \argmax{\pi(\theta|y,\tau)}
  fits[3] += -0.5*count_x*log(2*PI*Pmode[0]) - 0.5*sum_xx/Pmode[0];

  //maximised (log-)segment likelihood evaluated:
  fits[4] += -0.5*count_x*log(2*PI*sum_xx/count_x) - 0.5*count_x;

  free(Pmode);
  return;
}



//--------------------------------------------------------------------------
// Generic get functions
//--------------------------------------------------------------------------

void Get_Mprior(char **Mdist, Mprior_Ptr **Mprior, int *err){
  //Rprintf("Get_Mprior\n");
  if(      strcmp(*Mdist, "unif") == 0){
    *Mprior = Mprior_unif;
  }else if(strcmp(*Mdist, "pois") == 0){
    *Mprior = Mprior_pois;
  }else{
    *Mprior = Mprior_BLANK;
    *err = 1;
    return;
  }
}

void Get_Pprior(char **Pdist, Samp_Dist_Ptr **Samp_Dist, int *err){
  //Rprintf("Get_Pprior\n");

  if(      strcmp(*Pdist, "bern") == 0){
    *Samp_Dist = Samp_Dist_bern;
//  }else if(strcmp(*Pdist, "binom") == 0){
//    Samp_Dist = Samp_Dist_binom;
  }else if(strcmp(*Pdist, "pois") == 0){
    *Samp_Dist = Samp_Dist_pois;
  }else if(strcmp(*Pdist, "mean") == 0){
    *Samp_Dist = Samp_Dist_mean;
  }else if(strcmp(*Pdist, "norm") == 0){
    *Samp_Dist = Samp_Dist_norm;
  }else if(strcmp(*Pdist, "var" ) == 0){
    *Samp_Dist = Samp_Dist_var;
  }else{
    *Samp_Dist = Samp_Dist_BLANK;
    *err = 2;
    return;
  }
  return;
}

void Get_SumStats_FN(char **Pdist, Summary_Stats_Ptr **Summary_Stats, int *err){
  //Rprintf("Get_SumStats_FN\n");

  if(      strcmp(*Pdist, "bern") == 0){
    *Summary_Stats = Summary_Stats_bern;
//  }else if(strcmp(*Pdist, "binom") == 0){
//    Summary_Stats = Summary_Stats_binom;
  }else if(strcmp(*Pdist, "pois") == 0){
    *Summary_Stats = Summary_Stats_pois;
  }else if(strcmp(*Pdist, "mean") == 0){
    *Summary_Stats = Summary_Stats_mean;
  }else if(strcmp(*Pdist, "norm") == 0){
    *Summary_Stats = Summary_Stats_norm;
  }else if(strcmp(*Pdist, "var" ) == 0){
    *Summary_Stats = Summary_Stats_var;
  }else{
    *Summary_Stats = Summary_Stats_BLANK;
    *err = 3;
    return;
  }
  return;
}


void Get_SuffStats_Function(char **Pdist,
    Sufficient_Stats_Ptr **Sufficient_Stats, int *err){
  //Rprintf("Get_SuffStats_Function\n");

  if(      strcmp(*Pdist, "bern") == 0){
    *Sufficient_Stats = Sufficient_Stats_bern;
//  }else if(strcmp(*Pdist, "binom") == 0){
//    Sufficient_Stats = Sufficient_Stats_binom;
  }else if(strcmp(*Pdist, "pois") == 0){
    *Sufficient_Stats = Sufficient_Stats_pois;
  }else if(strcmp(*Pdist, "mean") == 0){
    *Sufficient_Stats = Sufficient_Stats_mean;
  }else if(strcmp(*Pdist, "norm") == 0){
    *Sufficient_Stats = Sufficient_Stats_norm;
  }else if(strcmp(*Pdist, "var" ) == 0){
    *Sufficient_Stats = Sufficient_Stats_var;
  }else{
    *Sufficient_Stats = Sufficient_Stats_BLANK;
    *err = 4;
    return;
  }
  return;
}

void Get_Param_Mode(char **Pdist, Param_Mode_Ptr **Param_Mode, int *err){
  //Rprintf("Get_Fit_FN\n");

  if(      strcmp(*Pdist, "bern") == 0){
    *Param_Mode = Param_Mode_bern;
//  }else if(strcmp(*Pdist, "binom") == 0){
//    Param_Mode = Param_Mode_binom;
  }else if(strcmp(*Pdist, "pois") == 0){
    *Param_Mode = Param_Mode_pois;
  }else if(strcmp(*Pdist, "mean") == 0){
    *Param_Mode = Param_Mode_mean;
  }else if(strcmp(*Pdist, "norm") == 0){
    *Param_Mode = Param_Mode_norm;
  }else if(strcmp(*Pdist, "var" ) == 0){
    *Param_Mode = Param_Mode_var;
  }else{
    *Param_Mode = Param_Mode_BLANK;
    *err = 5;
    return;
  }
  return;
}

void Get_Fit_FN(char **Pdist, Fit_FN_Ptr **Fit_FN, int *err){
  //Rprintf("Get_Fit_FN\n");

  if(      strcmp(*Pdist, "bern") == 0){
    *Fit_FN = Fit_FN_bern;
//  }else if(strcmp(*Pdist, "binom") == 0){
//    Fit_FN = Fit_FN_binom;
  }else if(strcmp(*Pdist, "pois") == 0){
    *Fit_FN = Fit_FN_pois;
  }else if(strcmp(*Pdist, "mean") == 0){
    *Fit_FN = Fit_FN_mean;
  }else if(strcmp(*Pdist, "norm") == 0){
    *Fit_FN = Fit_FN_norm;
  }else if(strcmp(*Pdist, "var" ) == 0){
    *Fit_FN = Fit_FN_var;
  }else{
    *Fit_FN = Fit_FN_BLANK;
    *err = 6;
    return;
  }
  return;
}



#endif //FILE_CASEFNS
