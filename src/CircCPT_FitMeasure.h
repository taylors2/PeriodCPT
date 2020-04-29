
/* 
//#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\stdio.h>
//#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\stdlib.h>
//#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\math.h>
#include <C:\Program Files\R\R-3.6.1\include\R.h>
//#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\time.h>
//#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\stdint.h>
//#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\string.h>
#include"C:\Users\staylo16\Desktop\CircCPT\src\CircCPT_General.h"
#include"C:\Users\staylo16\Desktop\CircCPT\src\CircCPT_MCMCcache.h"
#include"C:\Users\staylo16\Desktop\CircCPT\src\CircCPT_LOOKUP.h"
*/
 
 


int ModeCompareMCMCitem(MCMCitem_t *test, MCMCitem_t *instance){
  int fnname = 1005;
  profile(TRUE, fnname);
  int i;
  if(test->m != instance->m) goto DIFF;
  for(i=0;i<test->m;i++){
    if(test->tau[i] != instance->tau[i]) goto DIFF;
  }
  profile(FALSE, fnname);
  return TRUE;
  
  DIFF:;
  profile(FALSE, fnname);
  return TRUE;  //ERROR?!
}

MCMCitem_t * ModeTau(chain_t *chain1, chain_t *chain2, int *chainlen){
  int fnname = 1000;
  profile(TRUE, fnname);
  
  MCMCitem_t *draw;
  MCMCitem_t **address;
  int *tally;
  int nunique, i;

  address = (MCMCitem_t**)my_calloc(2 * *chainlen, sizeof(MCMCitem_t*));
  tally   = (int*)my_calloc(2 * *chainlen, sizeof(int));
  
  //initial from first chain
  draw = chain1->first;
  address[0] = draw;
  tally[0]++;
  draw = draw->next;
  nunique = 1;
  
  //run through first chain
  while(draw != NULL){
    for(i=0; i<nunique; i++){
      if(ModeCompareMCMCitem(draw, address[i])){
        tally[i]++;
        goto FOUNDSAME1;
      }
    }
    //draw is new
    address[nunique] = draw;
    tally[nunique]++;
    nunique++;
    
    
    FOUNDSAME1:;
    draw = draw->next;
  }
  
  //inital from second chain
  draw = chain2->first;
  
  //run through second chain
  while(draw != NULL){
    for(i=0; i<nunique; i++){
      if(ModeCompareMCMCitem(draw, address[i])){
        tally[i]++;
        goto FOUNDSAME2;
      }
    }
    //draw is new
    address[nunique] = draw;
    tally[nunique]++;
    nunique++;
    
    
    FOUNDSAME2:;
    draw = draw->next;
  }
  
  //find mode (returns first instance in the case of 2 equal modes!!!)
  int maxcount, modeplace;
  maxcount = 0;
  modeplace = 0;
  for(i=0; i<nunique; i++){
    if(tally[i]>maxcount){
      maxcount = tally[i];
      modeplace = i;
    }
  }
  
  MCMCitem_t *mode;
  mode = CopyMCMCitem(address[modeplace]);
  
  my_free(tally);
  my_free(address);
   
  profile(FALSE, fnname);
  return mode;
}



double **Normal_meanvar_SS(MCMCitem_t *TauHat, int *count, double *sum, double *sq, double *hypparam){
  
  int fnname = 1002;
  profile(TRUE, fnname);
  
  double **out;
  out = (double **)my_calloc(TauHat->m, sizeof(double*));
  int i;
  for(i=0;i<TauHat->m;i++){
    double *SS = (double *) my_calloc(4, sizeof(double));
    SS[1] = hypparam[1] + count[i]; //precision
    SS[0] = (hypparam[0]*hypparam[1] + sum[i]) / SS[1]; //mean
    SS[2] = hypparam[2] + 0.5*count[i]; //shape
    double xbar = sum[i]/count[i];
    double excess = xbar - hypparam[0];
    SS[3] = hypparam[3] + 0.5*(sq[i] - sum[i]*xbar) + 0.5*count[i]*hypparam[1]*excess*excess/SS[1]; //scale
    out[i] = SS;
  }
  
  profile(FALSE, fnname);
  return out;
  
}

void Normal_meanvar_Mode(MCMCitem_t *TauHat, double **SS, double *mode){
  
  int fnname = 1003;
  profile(TRUE, fnname);
  
  int i;
  for(i = 0; i<TauHat->m; i++){
    mode[i] = SS[i][0];
    mode[i + TauHat->m] = SS[i][3]/(SS[i][2] + 1.5);
      //Note: this calculates the ***variance*** mode rather than the precision!!!
  }
  
  profile(FALSE, fnname);
}


void Normal_meanvar_Fit(MCMCitem_t *TauHat, double **PostSS, int *count, double *sum, double *sq,
                           int *len, int *N, double *hypparam, double *hypS, double *hypM, int *minseglen,
                           int *maxM, double *ParamMode, double *fits){
  
  int fnname = 1004;
  profile(TRUE, fnname);
  
  int i;
  double tmp0;
  tmp0 = dTpois(TauHat->m, *maxM, *hypM, TRUE); //prior for changepoint number
  if(TauHat->m > 1){
    tmp0 -= log(*N);  //prior for 'a' changepoint event
    //prior for changepoint positions
    tmp0 += lgammafn(1 + *N - *minseglen * TauHat->m) + lgammafn(TauHat->m * *hypS) - lgammafn(*N - *minseglen * TauHat->m + *hypS * TauHat->m) - TauHat->m * lgammafn(*hypS);
    for(i=0; i<TauHat->m;i++){
      tmp0 += lgammafn(len[i] - *minseglen + *hypS) - lgammafn(1 + len[i] - *minseglen);
    }
  }
  
  //double MJTmode = tmp0;
  double JTmode  = tmp0;
  double MLmode  = 0;
  double Lmode   = 0;
  for(i=0; i<TauHat->m; i++){
    double muhat = ParamMode[i];
    double s2hat = ParamMode[i + TauHat->m];
    double mpost = PostSS[i][0];
    double ppost = PostSS[i][1];
    double apost = PostSS[i][2];
    double bpost = PostSS[i][3];

    //marginalised joint (like+prior) at TauHat.
    //MJTmode += -count[i]*0.5*log(2*PI) + hypparam[2]*log(hypparam[3]) + 0.5*log(hypparam[1]) - lgammafn(hypparam[2]) - apost*log(bpost) - 0.5*log(ppost) + lgammafn(apost);
    //joint (like+prior) at TauHat and mode estimate of segment parameters
    JTmode  += -(count[i]+1)*0.5*log(2*PI) + 0.5*log(hypparam[1]) + hypparam[2]*log(hypparam[3]) - 
      lgammafn(hypparam[2]) - (apost+1.5)*log(s2hat) - bpost/s2hat - 0.5*ppost*(muhat - mpost)*(muhat - mpost)/s2hat;
    //marginalised likelihood (wrt segment parameters) at TauHat
    MLmode  += -count[i]*0.5*log(2*PI) + hypparam[2]*log(hypparam[3]) + 0.5*log(hypparam[1]) - 
      lgammafn(hypparam[2]) - apost*log(bpost) - 0.5*log(ppost) + lgammafn(apost);
    //mode of likelihood at TauHat and mode estimate of segment parameters
    Lmode   += -count[i]*0.5*log(2*PI*s2hat) - 0.5*(sq[i] - 2*muhat*sum[i] + count[i]*muhat*muhat)/s2hat;
  }
    
  fits[0] = TauHat->value * -2;  // = MJTmode * -2;
  fits[1] = JTmode  * -2;
  fits[2] = MLmode  * -2;
  fits[3] = Lmode   * -2;
  
  profile(FALSE, fnname);
}



void Evaluate_Fit_and_SegParamModes(MCMCitem_t *TauHat, double *data, int *time,
                                    int *n, int *N, int *dist, double *hypparam,
                                    double *hypS, double *hypM, int *minseglen, int *Mmax,
                                    double *FitMeasure, double *ModeEst, int *nModeEst){
  
  int fnname = 1001;
  profile(TRUE, fnname);
  
  int *count, *len;
  int i;
  double *sum, *sq;
  int m = TauHat->m;
  
  len   = (int*)my_calloc(m, sizeof(int));
  count = (int*)my_calloc(m, sizeof(int));
  sum   = (double*)my_calloc(m, sizeof(double));
  sq    = (double*)my_calloc(m, sizeof(double));

  int iprev;
  for(i=0;i<m;i++){
    if(i==0){
      iprev = m - 1;
      len[i] = TauHat->tau[i] - TauHat->tau[iprev] + *N;
    }else{
      iprev = i - 1;
      len[i] = TauHat->tau[i] - TauHat->tau[iprev];
    }
  }
  
  int segID, j;
  for(i=0; i<*n; i++){
    segID = 0;
    for(j=0;j<m;j++){
      if(time[i] > TauHat->tau[j]) segID++;
    }
    if(segID == m) segID = 0;
    count[segID]++;
    sum[segID] += data[i];
    sq[segID] += data[i]*data[i];
  }  
  
  double **PostSS;
  switch(*dist){
    case 0:  //Poisson model -- hypparam = (alpha, beta)
      break;
    case 1:  //Binary model -- hypparam = (alpha, beta)
      break;
    case 2:  //Normal (mean prec) model -- hypparam = (m, c, a, b)
      PostSS = Normal_meanvar_SS(TauHat, count, sum, sq, hypparam);
      Normal_meanvar_Mode(TauHat, PostSS, ModeEst);
      Normal_meanvar_Fit(TauHat, PostSS, count, sum, sq, len, N, hypparam, hypS, hypM, minseglen, Mmax, ModeEst, FitMeasure);
      for(j=0; j<m; j++) my_free(PostSS[j]);
      my_free(PostSS);
      break;
    case 3:  //Mean model -- hypparam = (tau, m, c)
      break;
    case 4:  //Var model -- hypparam = (mu, a, b)
      break;
  }

  my_free(count);
  my_free(len);
  my_free(sum);
  my_free(sq);
  
  profile(FALSE, fnname);
}

void CopyModeFitToOutput(MCMCitem_t *TauHat, double *FitMeasure, double *ModeEst, int *nModeEst, int *out_mmode, int *out_taumode, double *out_segparammode, double *out_fit){
  int fnname = 1006;
  profile(TRUE, fnname);

  //copy to output

  *out_mmode = TauHat->m;
  int i;
  for(i = 0; i < TauHat->m; i++) out_taumode[i] = TauHat->tau[i];
  for(i = 0; i < *nModeEst; i++) out_segparammode[i] = ModeEst[i];
  for(i = 0; i < 4; i++) out_fit[i] = FitMeasure[i];

  //free memory
  my_free(ModeEst);
  my_free(FitMeasure);
  DeleteMCMCitem(TauHat);
  
  profile(FALSE, fnname);
}
  
