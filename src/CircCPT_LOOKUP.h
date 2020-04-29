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
*/

double Get_g2(int tau, int len, struct G2 *g2, int *N){
  //order: t1l1 t1l2 t1l3 t2l1 t2l2 t2l3 t3l1 t3l2 t3l3
  //Includes -Inf cases where len<minseglen (retain for spacing!!!)
  int fnname = 100;
  profile(TRUE, fnname);
  profile(FALSE, fnname);
  return g2[(tau-1) * *N + len - 1].value;
}

double Get_g1(int m, double *g1){
  int fnname = 101;
  profile(TRUE, fnname);
  profile(FALSE, fnname);
  return g1[m-1];
}


double dmhg(double *k, int length,  double K, double n, double N, int log_p){
  int fnname = 102;
  profile(TRUE, fnname);
  
  //Assume valid input values 
  double out = - lchoose(N,n);
  int i;
  for(i=0; i<length; i++)  out += lchoose(K,k[i]);
  if(log_p != TRUE)  out = exp(out);
  
  profile(FALSE, fnname);
  return out;
}


double dbbinom(double k, double n, double alpha, double beta, int log_p){
  int fnname = 103;
  profile(TRUE, fnname);
  
  //Assume valid input values
  double out;
  out = lgammafn(1 + n) - lgammafn(k+1) - lgammafn(n-k+1);
  out += lgammafn(alpha+beta) - lgammafn(alpha) - lgammafn(beta);
  out += lgammafn(k+alpha) + lgammafn(n-k+beta) - lgammafn(n+alpha+beta);
  //  double out = lchoose(n, k) + lbeta(k+alpha, n-k+beta) - lbeta(alpha, beta);
  if(log_p != TRUE) out = exp(out);
  
  profile(FALSE, fnname);
  return out;
}


double dmult( double *x, int length, double size, int log_p){
  int fnname = 104;
  profile(TRUE, fnname);

  //Assume valid input values
  int i;
  double probi = 1/((double) length);
  double out = lgammafn(size + 1);
  for(i=0; i<length; i++) out += x[i] * log(probi) - lgammafn(x[i]+1);
  if( log_p != TRUE ) out = exp(out);

  profile(FALSE, fnname);
  return out;
}

double dmvstMeanVar( double n, double x, double xx, double m, double c, double a, double b, int log_p){
  int fnname = 112;
  profile(TRUE, fnname);

  double a2 = a + n/2;
  double div = n+c;
  double mean = x/n;
  double cnst = lgammafn(a2) - lgammafn(a) - 0.5*n*log(2*PI) + 0.5*log(c) - 0.5*log(n+c) + a*log(b);
  double gr1 = 0.5*n*c*(mean-m)*(mean-m)/div;
  double gr2 = 0.5*(xx - 2*x*mean + n*mean*mean);
  double lgr = log(b + gr1 + gr2);
  double value = cnst - a2*lgr;
  if(log_p != TRUE) value = exp(value);

  profile(FALSE, fnname);
  return value;
}

double dmvstvar( double n, double x, double xx, double mu, double a, double b, int log_p){
  int fnname = 113;
  profile(TRUE, fnname);

  double a2 = a + n/2;
  double cnst = lgammafn(a2) - lgammafn(a) - 0.5*n*log(2*PI*a);
  double sc = 0.5*n*(log(a)-log(b));
  double gr2 = xx -2*x*mu + n*mu*mu;
  double lgr = log(1 + (gr2)/(2*b));
  double value = cnst + sc - a2*lgr;
  if(log_p != TRUE) value = exp(value);

  profile(FALSE, fnname);
  return value;
}

double dmvnorm( double n, double x, double xx, double tau, double m, double c, int log_p){
  int fnname = 114;
  profile(TRUE, fnname);

  double pt1 = -tau*0.5*(xx - x*x/n);
  double pt2 = -tau*n*c*0.5*(x*x/(n*n) - 2*m*x/n + m*m)/(n*tau+c);
  double cnst = 0.5*n*(log(tau)-log(2*PI)) + 0.5*(log(c) - log(n*tau+c));
  double value = cnst + pt1 + pt2;
  if(log_p != TRUE) value = exp(value);

  profile(FALSE, fnname);
  return value;
}


void ProcessData(double *data, int *time, int *n, int *N, double *sum, double *sq, int *count, double *orddata){
  int fnname = 105;
  profile(TRUE, fnname);
  
  int i;
  for(i=0; i<*n; i++){
    sum[ time[i] - 1] += data[i];
    sq[ time[i] - 1] += data[i]*data[i];
    count[ time[i] - 1]++;
  }
  for(i=0;i<*N;i++){ // duplicate
    sum[i+*N] = sum[i];
    sq[i+*N] = sq[i];
    count[i+*N] = count[i];
  }
  
  int cumcount[*N], skip[*N];
  cumcount[0] = 0;
  skip[0] = 0;
  for(i=1;i<*N;i++){
    cumcount[i] = cumcount[i-1] + count[i-1]; //#obs where 0<time<i
    skip[i] = 0; //#obs with time=(i+1) that has been ordered so far
  }

  int id;
  for(i=0;i<*n;i++){
    id = cumcount[time[i]-1] + skip[time[i]-1];
    orddata[id] = data[i];
    orddata[id+*n] = orddata[id]; //duplicate
    skip[time[i]-1]++;
  }

  profile(FALSE, fnname);
  return;
}



void Evaluate_G2(double *data, int *time, int *n, int *N, int *dist, double *hypparam,
                 double *hypS, int *minseglen, struct G2 *g2){
  int fnname = 106;
  profile(TRUE, fnname);

  //Pre-process data
  int *count;
  double *sq, *sum, *orddata;
  sum = (double *)my_calloc(2 * *N, sizeof(double));
  sq = (double *)my_calloc(2 * *N, sizeof(double));
  count = (int *)my_calloc(2 * *N, sizeof(int));
  orddata = (double *)my_calloc(2 * *n, sizeof(double));
  ProcessData(data, time, n, N, sum, sq, count, orddata);

  //Evaluate lookup table for function g2(tau, len)
  int len, tau, id, thisCount, startSeg, dataoffset, i;
  double thisSum, thisSq, val0, val1, val2;
  double Spart[*N];
  
  id = 0;
  for(tau=1; tau<=*N; tau++){
    thisSum = 0.0;
    thisSq = 0.0;
    thisCount = 0;

    //Find ordered data offset for time bin tau
    dataoffset = 0;
    i = 0;
    while(i < tau){
      dataoffset += count[i];
      i++;
    }

    for(len=1; len<=*N; len++){
      
      if(tau==1){ //calc & store the per segment contribution to the excess spacing between segments
        if((len<*N) & (len >= *minseglen)){
          Spart[len-1] = lgammafn((double)(len - *minseglen) + *hypS) ;
          Spart[len-1] -= lgammafn((double)(len - *minseglen) + 1);
          Spart[len-1] -= lgammafn( *hypS );
        }else{
          Spart[len-1] = 0.0;
        }
      }
      
      //Evaluate sum and count for segment (i.e. update last calc with extra obs at segStart)
      startSeg = tau - len + 1;
      if(startSeg<=0) startSeg += *N;   //range: 1--N
      thisSum += sum[startSeg - 1];
      thisSq += sq[startSeg - 1];
      thisCount += count[startSeg - 1];
      dataoffset -= count[startSeg - 1];  //move orddata offset back by the #obs in startSeg
      if(dataoffset<0) dataoffset += *n;
      
      
      if(len < *minseglen){  //Immediately deal with lengths shorter than minseglen
        g2[id].cpt = tau;
        g2[id].seglen = len;
        g2[id].value = -INFINITY; 
        id++;
        continue;
      }
      
      switch(*dist){
        case 0:  //Poisson model -- hypparam = (alpha, beta)
          val0 = Spart[len-1];
          val1 = dnbinom(thisSum, hypparam[0], hypparam[1]/(hypparam[1] + (double) thisCount), TRUE);
          val2 = dmult(orddata+dataoffset, thisCount, thisSum, TRUE);
          g2[id].cpt = tau;
          g2[id].seglen = len;
          g2[id].value = val0+val1+val2;
          break;
        case 1:  //Binary model -- hypparam = (alpha, beta)
          val0 = Spart[len-1];
          val1 = dbbinom(thisSum, (double) thisCount, hypparam[0], hypparam[1], TRUE);
          val2 = dmhg(orddata+dataoffset, thisCount, 1.0, thisSum, (double) thisCount, TRUE);
          g2[id].cpt = tau;
          g2[id].seglen = len;
          g2[id].value = val0 + val1 + val2;
          break;
        case 2:  //Normal (mean prec) model -- hypparam = (m, c, a, b)
          val0 = Spart[len-1];
          val1 = dmvstMeanVar((double)thisCount, thisSum, thisSq, hypparam[0], hypparam[1],
            hypparam[2], hypparam[3], TRUE);
          val2 = 0;
          g2[id].cpt = tau;
          g2[id].seglen = len;
          g2[id].value = val0 + val1 + val2;
          break;
        case 3:  //Mean model -- hypparam = (tau, m, c)
          val0 = Spart[len-1];
          val1 = dmvnorm((double)thisCount, thisSum, thisSq, hypparam[0], hypparam[1],
            hypparam[2], TRUE);
          val2 = 0;
          g2[id].cpt = tau;
          g2[id].seglen = len;
          g2[id].value = val0 + val1 + val2;
          break;
        case 4:  //Var model -- hypparam = (mu, a, b)
          val0 = Spart[len-1];
          val1 = dmvstvar((double)thisCount, thisSum, thisSq, hypparam[0], hypparam[1],
            hypparam[2], TRUE);
          val2 = 0;
          g2[id].cpt = tau;
          g2[id].seglen = len;
          g2[id].value = val0 + val1 + val2;
          break;
      }

      id++;
    }
  }

  my_free(sum);
  my_free(sq);
  my_free(count);
  my_free(orddata);
  
  profile(FALSE, fnname);
  return;
}


double dTpois(int m, int max, double rate, int log_p){
  int fnname = 107;
  profile(TRUE, fnname);
  
  //Truncated poisson(rate) density with range 1--max
  double out, area;
  area = ppois((double) max, rate, TRUE, FALSE) - ppois(0.0, rate, TRUE, FALSE);
  if(log_p == TRUE) area = log(area);
  out = dpois((double) m, rate, log_p);
  if(log_p == TRUE){
    out -= area;
  }else{
    out /= area;
  }
  
  profile(FALSE, fnname);
  return out;
}



void Evaluate_G1(int *maxM, int *N, double *hypS, double *hypM, int *minseglen, double *g1){
  int fnname = 108;
  profile(TRUE, fnname);
  
  int m = 1;
  g1[0] = dTpois(m, *maxM, *hypM, TRUE); //Only model prior for m=1, no changepoints!
  if(*maxM==1) {
    profile(FALSE, fnname);
    return;
  }
  
  for(m = 2; m <= *maxM; m++){
    //model prior
    g1[m-1] = dTpois(m, *maxM, *hypM, TRUE);	//Prior for number of segments
    g1[m-1] -= log((double) *N);  		//Prior for anchor cpt(?)
    //Global contribution to the excess spacing between segments
    g1[m-1] += lgammafn( (double) (*N - (m * *minseglen) + 1) );
    g1[m-1] += lgammafn( (double)m * *hypS);
    g1[m-1] -= lgammafn( (double)(*N - (m * *minseglen)) + ((double)m * *hypS));
  }
  
  profile(FALSE, fnname);
  return;
}  


void MAKE_LOOK_TABLES(double *data, int *time, int *n, int *N, int *dist, double *hypparam, double *hypS,
                      double *hypM, int *minseglen, int *maxM, double *g1, struct G2 *g2){
  int fnname = 109;
  profile(TRUE, fnname);  

  Evaluate_G1(maxM, N, hypS, hypM, minseglen, g1);
  Evaluate_G2(data, time, n, N, dist, hypparam, hypS, minseglen, g2);

  profile(FALSE, fnname);
  return;
}


double Eval_at_CPT(int *tau, int m, int *N, double *g1, struct G2 *g2){
  int fnname = 110;
  profile(TRUE, fnname);

  double out,tmp;
  int j;
  out = Get_g1(m,g1);
  if(m>1){
    tmp = Get_g2(tau[0], tau[0]-tau[m-1]+*N, g2, N);
    //tmp -= log((double)(tau[0]-tau[m-1]+*N));  //anchor first event based on d[0]. (WRONG)
    out += tmp;
    for(j = 1; j<m; j++){
      tmp = Get_g2(tau[j], tau[j]-tau[j-1], g2, N);
      out += tmp;
    }
  }else{
    tmp = Get_g2(tau[0], *N, g2, N);
    out += tmp;
  }
  profile(FALSE, fnname);
  return out;
}



void PrintLookup(double *g1, struct G2 *g2, int *maxM, int *N){
  int fnname = 111;
  profile(TRUE, fnname);


  FILE *f;
  f = fopen("g1table.txt","w");
  int id;
  fprintf(f,"model, g1\n");
  for(id=0; id<*maxM; id++){
    fprintf(f,"%d, %f\n",id+1,g1[id]);
  }
  fclose(f);
  f = fopen("g2table.txt","w");
  fprintf(f,"tau, seglen, value\n");
  for(id=0;id<(*N * *N);id++){
    fprintf(f,"%d, %d, %f\n", g2[id].cpt, g2[id].seglen, g2[id].value);
  }
  fclose(f);
  
  profile(FALSE, fnname);
  return;
}



