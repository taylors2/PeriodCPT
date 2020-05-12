#ifndef FILE_DISTRIBUTIONS
#define FILE_DISTRIBUTIONS

#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <math.h>


double dTpois(int m, int max, double rate, int log_p){
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
  return out;
}

double dmhg(double sumlchoose, double n, double N, int log_p){
  double out = sumlchoose - lchoose(N,n);
  if(log_p != TRUE)  out = exp(out);
  return out;
}

double dmult(double count, double sum, double sum_lfactorial, int log_p){
  //Assume prob vector is the same for all elements
  double out = lgammafn(sum + 1) - sum*log(count) - sum_lfactorial;
  if( log_p != TRUE ) out = exp(out);
  return out;
}

double dbbinom(double sum, double size, double alpha, double beta, int log_p){
  //Assume valid input values
  double out;
  out  = lgammafn(size + 1)            - lgammafn(sum + 1)     - lgammafn(size - sum + 1);
  out += lgammafn(alpha + beta)        - lgammafn(alpha)       - lgammafn(beta);
  out -= lgammafn(size + alpha + beta) - lgammafn(sum + alpha) - lgammafn(size - sum + beta);
  if(log_p != TRUE) out = exp(out);
  return out;
}

double dmvnorm( double n, double x, double xx, double precision, double loc, double sq_rate, int log_p){
  double pt1 = -precision*0.5*(xx - x*x/n);
  double pt2 = -precision*n*sq_rate*0.5*(x*x/(n*n) - 2*loc*x/n + loc*loc)/(n*precision+sq_rate);
  double cnst = 0.5*n*(log(precision)-log(2*PI)) + 0.5*(log(sq_rate) - log(n*precision+sq_rate));
  double value = cnst + pt1 + pt2;
  if(log_p != TRUE) value = exp(value);
  return value;
}

double dmvstMeanVar( double n, double x, double xx, double loc, double prec, double a, double b, int log_p){
  double a2 = a + n/2;
  double div = n + prec;
  double mean = x/n;
  double cnst = lgammafn(a2) - lgammafn(a) + a*log(b) + 0.5*(log(prec) - log(div) - n*log(2*PI));
  double gr1 = 0.5*n*prec*(mean-loc)*(mean-loc)/div;
  double gr2 = 0.5*(xx - n*mean*mean);
  double value = cnst - a2*log(b + gr1 + gr2);
  if(log_p != TRUE) value = exp(value);
  return value;
}

double dmvstvar( double n, double x, double xx, double mu, double a, double b, int log_p){
  double a2 = a + 0.5*n;
  double cnst = lgammafn(a2) - lgammafn(a) - 0.5*n*log(2*PI) + a*log(b);
  double value = cnst - a2*log(b + 0.5*xx - x*mu + 0.5*n*mu*mu);
  if(log_p != TRUE) value = exp(value);
  return value;
}






#endif //FILE_DISTRIBUTIONS
