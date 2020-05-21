\name{table_pcpt}
\alias{table_pcpt}

\title{
  ~~Frequency table of the within period changepoint events.~~
}
\description{
  Frequency table of the within period changepoint events.
}
\usage{
  table_pcpt(object)
}

\arguments{
  \item{object}{
    An pcpt class object.
  }
}

\details{
  Returns a frequency table of the number of samples where the within period time step appears anywhere in the periodic changepoint vector. Note that the changepoint vector "c(1)" corresponds to the no-within period change (m=1) where the "1" is used as a holding space. As such, these samples are added as an extra column at the end of the table and the totals elsewhere in the table relates to the remaining samples where the number of within period changepoint events are >= 2.


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
  table_pcpt(ans)
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
\keyword{changepoint}
