\name{table_npcpt}
\alias{table_npcpt}

\title{
  ~~Frequency table of the number of within period changepoint events.~~
}
\description{
  Frequency table of the number of within period changepoint events
}
\usage{
  table_npcpt(object)
}

\arguments{
  \item{object}{
    An pcpt class object.
  }
}

\details{
  Returns a frequency table of the number of within period changepoint events in the MCMC samples.
}

\value{
  A table either with each row corresponding to each chain, or a single row if the chains have previously been summarised into a single overall table.
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
  ans  = PeriodCPT.bern(data)
  ans  = summarise_chains(ans)
  table_npcpt(ans)
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
\keyword{changepoint}
