\name{Param_quantiles}
\alias{Param_quantiles}
\title{Estimate segment parameters}
\description{
  Estimate segment parameters per time index.
}
\usage{
  Param_quantiles(probs=0.5, DRAW)
}

\arguments{
  \item{DRAW}{Output from CircCPT.}
  \item{probs}{Vector of probabilities.}
}

\details{
  Evaluates the quantile estimates for the parameter (rate for 
Poisson model or probability for Binomial model) per time bin.

}

\value{
  A vector (or matrix if multiple probabilities are specified) 
for the quantile estimate for each time bin. If multiple 
probabilities are specified, then a matrix is returned instead 
with each column containing the vector quantile estimate for 
each time bin for a given probability.
}

\examples{
  set.seed(10)
  rate = c(rep(1,12),rep(10,12),rep(1,12),rep(10,12))
  data = rpois(48,rate)
  time = c(1:24,1:24)
  DRAW <- CircCPT(data=data, time=time)
  Q <- Param_quantiles(probs=c(0.5,0.025,0.975), DRAW)
  head(Q)
}

\keyword{Param_quantiles}
