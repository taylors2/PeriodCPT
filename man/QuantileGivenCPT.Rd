\name{QuantileGivenCPT}
\alias{QuantileGivenCPT}
\title{Estimate segment parameters}
\description{
  Evaluates the quantile estimates of the segment parameters for
a given vector of changepoint positions.
}
\usage{
  QuantileGivenCPT(probs=0.5, cpts = NULL, data = NULL,
    time = NULL, l = 4, hyp=NULL,
    DIST = c("Poisson","Binomial","Normal","Mean","Var"))
}

\arguments{
  \item{probs}{Vector of probabilities.}
  \item{cpts}{Vector of changepoint positions.}
  \item{data}{Data vector.}
  \item{time}{time vector.}
  \item{l}{minimum segment length.}
  \item{hyp}{Prior hyper parameters.}
  \item{DIST}{Distributional assumption.}
}

\details{
  Evaluates the quantile probabilities for the segment 
parameters for a given segmentation.

  See manual for "CircCPT" for full details of the 
distibution options.  
}

\value{
  Returns a vector (matrix if more than 1 probabilities are 
specified) of the quantile estimates for the parameter on each 
segment as specified by the given changepoint vector.
}


\examples{
  set.seed(10)
  rate = c(rep(1,12),rep(10,12),rep(1,12),rep(10,12))
  data = rpois(48,rate)
  time = c(1:24,1:24)
  QuantileGivenCPT(probs=c(0.5,0.025,0.975), cpts=c(12,24),
    data=data, time=time, DIST="Poisson")
}

\keyword{QuantileGivenCPT}
