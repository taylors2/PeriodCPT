///*
#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

//*/

/*
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\stdio.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\stdlib.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\math.h>
#include <C:\Program Files\R\R-3.4.3\include\R.h>
#include <C:\Program Files\R\R-3.4.3\include\Rmath.h>
#include <C:\Program Files\R\R-3.4.3\include\Rinternals.h>
#include <C:\Program Files\R\R-3.4.3\include\R_ext/Rdynload.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\time.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\stdint.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\string.h>

 */

extern void CircCPT_RJMCMC(
    double *data,       //data vector
    int *time,          //time vector
    int *n,             //length of data vector
    int *N,             //Period duration
    int *dist,          //Distribution assumption
    double *hypparam,   //parameter prior hyper parameters
    double *hypS,       //cpt spread hyper parameter
    double *hypM,       //Model prior hyper parameter
    int *minseglen,     //minimum segment length
    int *itstep,        //Number of MCMC iterations per batch
    int *Maxit,         //Number of MCMC batches
    double *typeI,      //Type I error rate for convergence test
    int *init1,         //initial cpt vector for chain 1
    int *m1,            //length of init1
    int *init2,         //initial cpt vector for chain 2
    int *m2,            //length of init2
    int *progress,      //print progress bar
    int *CACHEMAX,      //size of cache
    int *maxM,          //maximum length of changepoint vector
    int *burn,          //immediately ignore first batch as burn-in
    int *err,           //error flag (0=ok, 1=dist error, 2=ndraw too small)
    int *draw,          //output chain samples
    int *ndraw,         //size of output vecotor
    int *conv,           //Convergence satisfied
    int *numParamPerSeg,  //Number of segment parameters per segment
    int *out_mmode,             //Mode estimate for the number of changepoints
    int *out_taumode,           //Mode estimate for the changepoint positions
    double *out_segparammode,   //Mode estimates for the segment parameters
    double *out_fits            //Measures of fit
);


#include"CircCPT_General.h"
#include"CircCPT_LOOKUP.h"
#include"CircCPT_MCMCcache.h"
#include"CircCPT_PROPfunctions.h"
#include"CircCPT_MCMCfunctions.h"
#include"CircCPT_FitMeasure.h"

void CircCPT_RJMCMC(double *data, int *time, int *n, int *N, int *dist, double *hypparam,
                    double *hypS, double *hypM, int *minseglen, int *itstep, int *Maxit, double *typeI,
                    int *init1, int *m1, int *init2, int *m2, int *progress, int *CACHEMAX, int *maxM, int *burn,
                    int *err, int *draw, int *ndraw, int *conv, int *numParamPerSeg,
                    int *out_mmode, int *out_taumode, double *out_segparammode, double *out_fits){
  int fname = 000;
  profile(TRUE, fname);

  //Dist: 0=Poisson, 1=Binomial, 2=Normal(meanvar), 3=Norm(mean), 4=Norm(var)
  //time: in range 1,...,N
  //err: 1=invalid dist, 2=insufficient draw length to export all samples.
  //conv: 0=Failed to converge, 1=convergence successful.

  if((*dist < 0) | (*dist>4)){  //Is dist specified correctly
    *err = 1;
    goto ESCAPE;
  }

  //MAKE LOOKUP TABLES
  double *g1;
  struct G2 *g2;
  g1 = (double *) my_calloc( *maxM ,  sizeof(double) );
  g2 = (struct G2 *) my_calloc( *N * *N , sizeof(struct G2) );
  MAKE_LOOK_TABLES(data, time, n, N, dist, hypparam, hypS, hypM, minseglen, maxM, g1, g2);

  //Make & Initial MCMC chains by pushing inits
  chain_t *MCMC1, *MCMC2;
  MCMC1 = MakeMCMCchain();
  PushToChain(MCMC1, MakeMCMCitem(*m1, init1, N, g1, g2, 0));
  MCMC2 = MakeMCMCchain();
  PushToChain(MCMC2, MakeMCMCitem(*m2, init2, N, g1, g2, 0));


  //Perform batch of MCMC samples
  int iter = 0;
  int passTest = FALSE;
  int BURN = *burn;
  while(((*Maxit==0) | (iter < *Maxit))  & !passTest){
    iter++;

    //<<<\/>>> THIS CAN BE DONE IN PARALLEL <<<\/>>>

    //run MCMC 1
    ClearChain(MCMC1);
    batch_iteration(MCMC1, itstep, g1, g2, minseglen, N, progress, CACHEMAX);

    //run MCMC 2
    ClearChain(MCMC2);
    batch_iteration(MCMC2, itstep, g1, g2, minseglen, N, progress, CACHEMAX);

    //<<</\>>> THIS CAN BE DONE IN PARALLEL <<</\>>>

    if(BURN==FALSE){
      //preform convergence test
      passTest = ConvTest(MCMC1, MCMC2, typeI, N, maxM);
    }else{
      BURN = FALSE;
    }
  }

  *conv = passTest; // 0=non-convergence, 1=converge

  //***new: Evaluate mode estimats and calculate fit metrics
  int TotalChainLength = MCMC1->length + MCMC2->length;
  MCMCitem_t *TauHat = ModeTau(MCMC1, MCMC2, &TotalChainLength);
  double *FitMeasure = (double*)my_calloc(4, sizeof(double));
  int nModeEst = TauHat->m * *numParamPerSeg;
  double *ModeEst = (double*)my_calloc(nModeEst * 2, sizeof(double));
  Evaluate_Fit_and_SegParamModes(TauHat, data, time, n, N, dist, hypparam, hypS, hypM, minseglen, maxM,
                                 FitMeasure, ModeEst, &nModeEst);
  CopyModeFitToOutput(TauHat, FitMeasure, ModeEst, &nModeEst, out_mmode, out_taumode, out_segparammode, out_fits);
  //***new: end


  //export MCMC samples
  ExportDraws(MCMC1, MCMC2, draw, ndraw, err);

  //free allocated memory
  DeleteMCMCchain(MCMC2);
  DeleteMCMCchain(MCMC1);

  my_free(g2);  //free g2 lookup table
  my_free(g1);  //free g1 lookup table
  ESCAPE:;

  profile(FALSE, fname);
  return;
}


void R_init_CirCPT_Cfunctions(DllInfo *info){

  R_CMethodDef cMethods[] = {
    {"CircCPT_RJMCMC", (DL_FUNC) &CircCPT_RJMCMC, 29},
    {NULL,NULL,0}
  };

  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

