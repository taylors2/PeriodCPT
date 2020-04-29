% File src/library/base/man/CircCPT.Rd
\name{CircCPT}
\alias{CircCPT}
\title{Sample from changepoint posterior}
\description{
  Performs an RJMCMC sampler to approximate the periodic 
changepoint posterior distribution.
}
\usage{
  CircCPT(data, time, DIST=c("Poisson","Binomial","Normal",
    "Mean","Var"), hyp=NULL, gamma=1, delta=1, l=4, 
    init1=NULL, init2=NULL, control=list(), summary=TRUE,
    MODE = FALSE)
}

\arguments{
  \item{data}{Data vector.}
  \item{time}{Vector of integers with time bin index.}
  \item{DIST}{Either "Poisson" for count data or "Binomial" for 
binary data.}
  \item{hyp, gamma, delta}{Prior hyper parameters}
  \item{l}{Minimum segment length.}
  \item{init1, init2}{Integer vector giving the initial 
changepoint vector for the two MCMC chains. If NULL, then the 
chain is initiated from opposite side of paramter space.}
  \item{control}{A list containing the algorithm parameters, 
see details below.}
  \item{summary}{Logical flag. If TRUE then the draws are 
summarised. Otherwise the the chains are returned in order.}
  \item{MODE}{Logical flag. If TRUE then the MAP estimates 
  are evaluated as well as metric of fit. Default, FALSE.
  (Currently only working for when DIST="Normal")}
}

\details{
  Details regarding the different distribution assumptions:
  \itemize{
    \item Poisson: change in rate parameter in Poisson count 
      data. Prior distribution for the parameter is Gamma. 
      Argument 'hyp' contains the two prior hyperparameters, 
      specified as hyp=c(1,1) by default.
    \item Binomial: change in success probability parameter
     assuming a fixed size of 1. Prior distribution for the 
     parameter is Beta. Argument 'hyp' contains the two prior 
     hyperparameters, specified as hyp=c(1,1) by default.
    \item Normal: change in mean and variance in real valued
      data. Joint Normal-InvGamma prior distribution is defined
      for the mean & variance. Argument 'hyp' contains the 
      four prior hyperparameters (m,c,a,b), specified as 
      hyp=c(0,1,1,1) by default.
    \item Mean: change in mean in real valued data with assumed
      known constant variance. Prior distribution
      for the mean is Normal. Argument 'hyp' contains the known 
      variance followed by the Normal location (m) and scale (c)
      hyperparameters, sepcified as hyp=c(1,1,1) by default.
    \item Var: change in variance in real valued data with 
      assumed known constant mean. The distribution for the 
      variance is InvGamma. Argument 'hyp' contains the known 
      mean followed by the shape (a) and scale (b)
      hyperparameters, specified as hyp=c(0,1,1) by default.
  }
  Note: Variance in Normal, Mean and Variance assumptions are 
  transformed into precision for calculation ease and 
  subsequently reverted back on output.
  
  The control argument is a list with the following items:
  \itemize{
    \item itstep: number of iterations to perform in each batch 
(default: 1000)
    \item Maxit: Maximum number of batches to perform 
(default: 10). If NULL then sampler continues until convergence 
between the two chains is satisfied. Only the samples from the 
last batch is returned.
    \item cachesize: Maximum size of the cache in storing the 
proposal sets.(default: 50)
    \item typeI: Type I error when performing the convergence 
test. (default: 0.05)
    \item burn: Logical, should the first batch be discarded a 
burn-in. (default: TRUE)
    \item progress: Logical, should a progress bar per chain 
and batch be printed to screen (default: FALSE).
  }

}

\value{
  Returned is a list with the items TAU (a matrix containing 
the draws), data and time.
  
  If Summary=TRUE, then TAU is summarised such that each row 
contains a unique sample and the columns in order defined the 
number of segments (M), the frequency of the cpt vector (Freq), 
the changepoint vector (tau1,...) and the posterior statistics 
for the segment parameter (e.g. alpha1,... and beta1,... for 
Poisson distribution assumption).
  
  If Summary=FALSE, then TAU only contains the changepoint 
vectors with each row denoting the iteration index.
  
  The attributes of the returned list contain the statistical 
and algorithmic parameters as well as performance information.
  
  If MODE=TRUE, then an additional item named "MODE" is added to the
output which is a list containing the MAP estimates and fit metrics.
}

\examples{
  set.seed(10)
  rate = c(rep(1,12),rep(10,12),rep(1,12),rep(10,12))
  data = rpois(48,rate)
  time = c(1:24,1:24)
  DRAW <- CircCPT(data=data, time=time)
  attributes(DRAW)$Message
  DRAW$TAU
}

\keyword{CircCPT}
