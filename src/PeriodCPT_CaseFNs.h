#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <math.h>

#include "PeriodCPT_distributions.h"

//--General structure for functions
//double Mprior_"Mdist"(int m, int *max, double *params);
//  m      - pcpt count to evalute prior.
//  *max   - upper bound of m.
//  *param - distribution specific parameters.
//  RETURN - double of log-pmf of pcpt count prior evaluated at m.

//double **Summary_Stats_"Pdist"(double *data, int *time, int *n, int *N);
//  *data  - time series data
//  *time  - within period time instance corresponding to data
//  *n     - length of time series
//  *N     - period length
//  RETURN - pointer to pointer array of summary statistics, with first
//             element indicating the number of statistics in array.

//double Samp_Dist_"Pdist"(double *SumStats, int nStat, double *Phyp);
//  *SumStats - Summary statistics
//  nStat     - number of summary statistics, length of SumStats
//  *Phyp     - Prior hyper-params for segment parameters
//  RETURN    - log-marginal density/mass between sampling distribution and segment parameter priors.

//-------------------------------------------------------------

void Free_Summary_Stats(double **SumStats){
  double number_of_stat_types = SumStats[0][0];
  for(int i = 1; i < number_of_stat_types; i++){
    free(SumStats[i]);
  }
  free(SumStats[0]);
  free(SumStats);
  return;
}

double **Make_Summary_Stats(int number_of_stat_types, int length){
  double **SumStats;
  SumStats = (double **)calloc(1 + number_of_stat_types, sizeof(double *));
  SumStats[0] = (double *)calloc(1, sizeof(double));
  SumStats[0][0] = (double)number_of_stat_types;
  for(int i = 1; i <= number_of_stat_types; i++){
    SumStats[0] = (double *)calloc(2 * length, sizeof(double));
  }
  return SumStats;
}

void Get_Functions(char **Mdist, char **Pdist, double (*Mprior)(), double (*Samp_Dist)(),
                   double **(*Summary_Stats)(), int *err){

//double (*Mprior)();
  double Mprior_unif();
  double Mprior_pois();
  if(      strcmp(*Mdist, "unif") == 0){
    Mprior = Mprior_unif;
  }else if(strcmp(*Mdist, "pois") == 0){
    Mprior = &Mprior_unif;
  }else{
    *err = 1;
    return;
  }

//double   (*Samp_Dist)();        //returns a double
//double **(*Summary_Stats)();    //returns a double pointer containg double values
  double   Samp_Dist_bern();
  double **Summary_Stats_bern();
  double   Samp_Dist_binom();
  double **Summary_Stats_binom();
  double   Samp_Dist_mean();
  double **Summary_Stats_mean();
  double   Samp_Dist_norm();
  double **Summary_Stats_norm();
  double   Samp_Dist_pois();
  double **Summary_Stats_pois();
  double   Samp_Dist_var();
  double **Summary_Stats_var();

  if(      strcmp(*Pdist, "bern") == 0){
    Samp_Dist = &Samp_Dist_bern;
    Summary_Stats = &Summary_Stats_bern;
//  }else if(strcmp(*Pdist, "binom") == 0){
//    Samp_Dist = &Samp_Dist_binom;
//    Summary_Stats = &Summary_Stats_binom;
  }else if(strcmp(*Pdist, "pois") == 0){
    Samp_Dist = &Samp_Dist_pois;
    Summary_Stats = &Summary_Stats_pois;
  }else if(strcmp(*Pdist, "mean") == 0){
    Samp_Dist = &Samp_Dist_mean;
    Summary_Stats = &Summary_Stats_mean;
  }else if(strcmp(*Pdist, "norm") == 0){
    Samp_Dist = &Samp_Dist_norm;
    Summary_Stats = &Summary_Stats_norm;
  }else if(strcmp(*Pdist, "var" ) == 0){
    Samp_Dist = &Samp_Dist_var;
    Summary_Stats = &Summary_Stats_var;
  }else{
    *err = 2;
    return;
  }
  return;
}

//------------------------------------------------------------

double Mprior_pois(int m, int *max, double *params){
  //Trunc-Pois on interval 1, ..., max.
  double rate = params[0];
  return dTpois(m, *max, rate, TRUE);
}

double Mprior_unif(int m, int *max, double *params){
  //Discrete-Unif on interval 1, ..., max.
  return -log(*max);
}

//------------------------------------------------------------

double Samp_Dist_bern(double *SumStats, int nStat, double *Phyp){
  double sum_x          = SumStats[1];
  double count_x        = SumStats[2];
  double alpha          = Phyp[0];
  double beta           = Phyp[1];
  return dbbinom(sum_x, count_x, alpha, beta, TRUE) +
           dmhg(1.0, sum_x, count_x, TRUE);
}

double **Summary_Stats_bern(double *data, int *time, int *n, int *N){
  double **SumStats = Make_Summary_Stats(2, *N);
  for(int i = 0; i < *n; i++){
    SumStats[1][time[i] - 1] += data[i];    //1: sum_xi
    SumStats[2][time[i] - 1] += 1;          //2: count_xi
  }
  return SumStats;
};

//------------------------------------------------------------

/*
double Samp_Dist_binom(double *SumStats, int nStat, double *Phyp){
  double sum_x          = SumStats[1];
  double sum_n          = SumStats[2];
  double sum_lchoose_nx = SumStats[3];
  double alpha          = Phyp[0];
  double beta           = Phyp[1];
  return dbbinom(sum_x, sum_n, alpha, beta, TRUE) +
    dmhg(sum_lchoose_nx, sum_x, sum_n, TRUE);
}


double **Summary_Stats_binom(double *data, int *time, int *n, int *N){
  double **SumStats = Make_Summary_Stats(3, *N);
  for(int i = 0; i < *n; i++){
    SumStats[1][time[i] - 1] += data[i];                        //1: sum_xi
    SumStats[2][time[i] - 1] += data[i + *n];                   //2: sum_ni
    SumStats[3][time[i] - 1] += lchoose(data[i + *n], data[i]); //3: sum_log(choose(ni, xi))
  }
  return SumStats;
};
*/

//------------------------------------------------------------

double Samp_Dist_pois(double *SumStats, int nStat, double *Phyp){
  double sum_x          = SumStats[1];
  double count_x        = SumStats[2];
  double sum_lfact      = SumStats[3];
  double alpha          = Phyp[0];
  double beta           = Phyp[1];
  return dnbinom(sum_x, alpha, beta/(beta + count_x), TRUE) +
    dmult(count_x, sum_x, sum_lfact, TRUE);
}

double **Summary_Stats_pois(double *data, int *time, int *n, int *N){
  double **SumStats = Make_Summary_Stats(3, *N);
  for(int i = 0; i < *n; i++){
    SumStats[1][time[i] - 1] += data[i];            //1: sum_xi
    SumStats[2][time[i] - 1] += 1;                  //2: count_xi
    SumStats[3][time[i] - 1] += lgamma(data[i]+1);  //3: sum_log(xi!)
  }
  return SumStats;
};

//------------------------------------------------------------

double Samp_Dist_mean(double *SumStats, int nStat, double *Phyp){
  double count_x        = SumStats[1];
  double sum_x          = SumStats[2];
  double sum_xx         = SumStats[3];
  double m              = Phyp[0];
  double c              = Phyp[1];
  double var            = Phyp[2];
  return dmvnorm(count_x, sum_x, sum_xx, m, 1/c, 1/var, TRUE);
}

double **Summary_Stats_mean(double *data, int *time, int *n, int *N){
  double **SumStats = Make_Summary_Stats(3, *N);
  for(int i = 0; i < *n; i++){
    SumStats[1][time[i] - 1] += 1;                  //1: count_xi
    SumStats[2][time[i] - 1] += data[i];            //2: sum_xi
    SumStats[3][time[i] - 1] += data[i]*data[i];    //3: sum_xi^2
  }
  return SumStats;
};

//------------------------------------------------------------

double Samp_Dist_norm(double *SumStats, int nStat, double *Phyp){
  double count_x = SumStats[1];
  double sum_x   = SumStats[2];
  double sum_xx  = SumStats[3];
  double m       = Phyp[0];
  double c       = Phyp[1];
  double alpha   = Phyp[2];
  double beta    = Phyp[3];
  return dmvstMeanVar(count_x, sum_x, sum_xx, m, 1/c, alpha, beta, TRUE);
}

double **Summary_Stats_norm(double *data, int *time, int *n, int *N){
  double **SumStats = Make_Summary_Stats(3, *N);
  for(int i = 0; i < *n; i++){
    SumStats[1][time[i] - 1] += 1;                  //1: count_xi
    SumStats[2][time[i] - 1] += data[i];            //2: sum_xi
    SumStats[3][time[i] - 1] += data[i]*data[i];    //3: sum_xi^2
  }
  return SumStats;
};

//------------------------------------------------------------

double Samp_Dist_var(double *SumStats, int nStat, double *Phyp){
  double count_x = SumStats[1];
  double sum_xx  = SumStats[2];
  double alpha   = Phyp[2];
  double beta    = Phyp[3];
  return dmvstvar(count_x, 0.0, sum_xx, 0.0, alpha, beta, TRUE);
}

double **Summary_Stats_var(double *data, int *time, int *n, int *N){
  double **SumStats = Make_Summary_Stats(3, *N);
  for(int i = 0; i < *n; i++){
    SumStats[1][time[i] - 1] += 1;                  //1: count_xi
    SumStats[2][time[i] - 1] += data[i]*data[i];    //2: sum_xi^2
  }
  return SumStats;
};


