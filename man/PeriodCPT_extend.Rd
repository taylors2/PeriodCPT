\name{PeriodCPT_extend}
\alias{PeriodCPT_extend}

\title{
  ~~Run more MCMC samples.~~
}
\description{
  Run more MCMC samples.
}
\usage{
  PeriodCPT_extend(object, newiters = 1e4)
}

\arguments{
  \item{object}{
    An pcpt class object.
  }
  \item{newiters}{
    A positive integer indicating how many iterations to perform.
  }
}

\details{
This function provides the utility to run the with period cpt algorithm for a further \code{newiters} iterations following on from the last sample from the previous run. All previous samples will be treated as being part of the burn-in period and will be discarded.
}

\value{
  An object of class \code{pcpt} with the new chains.
}

\author{
Simon Taylor
}

\seealso{
  \code{\linkS4class{pcpt}}, \code{\link{PeriodCPT}}
}

\examples{
  set.seed(1)
  data = ts( rbinom(120, size = 1, prob = rep(c(0.2, 0.8), each=6)), freq = 12)
  ans  = PeriodCPT.bern(data, niter = 1e4)
  ans  = PeriodCPT_extend(ans, newiters = 1e4)
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
\keyword{changepoint}
