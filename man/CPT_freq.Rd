\name{CPT_freq}
\alias{CPT_freq}
\alias{CPT_mode}
\title{Marginal cpt position table}
\description{
  Estimates the probability mass function (pmf) for the 
marginal changepoint positions or joint changepoint vectors.
}
\usage{
  CPT_freq(DRAW, marg=TRUE)
  CPT_mode(DRAW)
}

\arguments{
  \item{DRAW}{Output from CircCPT.}
  \item{marg}{Logical flag, if TRUE then the marginal pmf is 
evaluated for the changpoint positions. If FALSE, then the 
joint pmf is evaluated for the changepoint vector.}
}

\details{
  Estimates the probability mass function (pmf) for the 
marginal changepoint positions or joint changepoint vectors.
}

\value{
  CPT_freq() returns a table whereas CPT_mode() returns a 
vector containing the times of the most probable changepoint 
vector.
}

\examples{
  set.seed(10)
  rate = c(rep(1,12),rep(10,12),rep(1,12),rep(10,12))
  data = rpois(48,rate)
  time = c(1:24,1:24)
  DRAW <- CircCPT(data=data, time=time)
  head(CPT_freq(DRAW))
  CPT_mode(DRAW)
}

\keyword{CPT_freq}
