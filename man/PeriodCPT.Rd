\name{PeriodCPT}
\alias{PeriodCPT}
\title{
  ~~function to evaluate within period changepoint events ~~
}
\description{
  Wrapper function for evaluating within period changepoint events for a specified sampling distribution.
}
\usage{
PeriodCPT(data, distribution = c("bern","pois","norm","norm","mean","var"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    A vector or ts object containing the data within which you which to find within period changepoint events. If not specified, the period length is taken from the ts object.
  }
  \item{distribution}{
    Choice of "bern", "pois", "norm", "mean" and "var" that specifies the sampling distribution and the type of statistic to investigate within period changepoint events.
  }
  \item{...}{
    Additional optional arguments that specifies prior details, MCMC parameters and output format.
  }

}

\details{
  This function is used to find with period changepoint events for selected sampling distribution and statistic. A changepoint is denoted as the (!)LAST OBSERVATION OF THE CURRENT SEGMENT/REGIME(!)

  This is a wrapper function for different distribition options. See the appropriate PeriodCPT.`distribution` for more information on implementation.
}

\value{
  An object of S4 class "pcpt" is returned.
}

\author{
Simon Taylor
}

\seealso{
  \code{\linkS4class{pcpt}}, \code{\link{PeriodCPT.bern}}, \code{\link{PeriodCPT.pois}}, \code{\link{PeriodCPT.norm}}, \code{\link{PeriodCPT.mean}}, \code{\link{PeriodCPT.var}}
}

\examples{
  data = ts( rbinom(120, size = 1, prob = rep(c(0.2, 0.8), each=6)), freq = 12)
  ans  = PeriodCPT(data, distribution = "bern")
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
\keyword{changepoint}
