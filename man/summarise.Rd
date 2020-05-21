\name{summarise}
\alias{summarise}
\alias{summarize}
\alias{summarise_chains}
\alias{summarize_chains}
\alias{summarise_combine}
\alias{summarize_combine}

\title{
  ~~functions to summarise the pcpt chain results.~~
}
\description{
  Procedures to summarise the MCMC results into a summary table(s) with sufficienct statistics.
}
\usage{
  summarise_chains(object, all = TRUE)
  summarize_chains(object, all = TRUE)
  summarise_combine(object)
  summarize_combine(object)
}

\arguments{
  \item{object}{
    An pcpt class object.
  }
  \item{all}{
    Logical indicating whether all of the chains should be combined into a single summary table.
  }
}

\details{
  The output of \code{PeriodCPT()} contains the full MCMC samples per chain within the results slot. The function \code{summarise_chains()} creates a frequency table where each row describes a unique sample. The order of the columns are:
  \describe{
    \item{m}{Number of within period changepoint events.}
    \item{freq}{A tally of the number of times the within period changepoint vector appears in the chain.}
    \item{tau[j]}{The within period changepoint event. An NA cell is used to pad out short changepoint vectors such that the the array can account for changepoint vector that have more events.}
  }
  The following columns contain the posterior sufficient statistics for the segment parameters for the the corresponding segments depending on the specified sampling distribution.

  If \code{all=FALSE}, then the table is evaluated for each chain, and the results slot will contain a list of tables in the appropriate order. Otherwise, if \code{all=TRUE} then the tables from each chain is combined into a single overall table. In this case, the results slot will contain a list with only one item of the overall table. The \code{summarise_combine()} function provides the utility to subsequently combine multiple summarised tables if \code{all = FALSE} when \code{summaise_chains()} was first called.

  Note that \code{summarise_chains()} changes logical flag of the summarise slot of the pcpt object to indicate that the results slot contains the summarised values rather than the raw chains. Also, the MCMC.last slot is populated with the last sample per chain to permit more MCMC samples to be drawn.

  The functions \code{summarize_chains()} and \code{summarize_combine()} implements the corresponding procedures.

}

\value{
  An object of S4 class "pcpt" is returned.
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
  head(results(ans)[[1]])
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
\keyword{changepoint}
