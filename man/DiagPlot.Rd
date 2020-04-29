\name{diagnostic.plot}
\alias{diagnostic.plot}
\title{
  diagnostic plot of posterior samples
}
\description{
  Diagnostic plot of the posterior samples obtained by 
CircCPT() function.
}

\usage{
  diagnostic.plot(DRAW, night=NULL, NAME)
}

\arguments{
  \item{DRAW}{Output from CircCPT.}
  \item{night}{Vector containing two items denoting the start 
and end of night times. If NULL, then night is ignored.}
  \item{NAME}{Plot title.}
}

\details{
  An image is produced with three frames:
  
  Top-left: most probable changepoint positions for a given 
number of segments.
  
  Bottom-right: barplot of the marginal posterior for the 
number of segments.
  
  Right: Circular plot presenting the data as well as the 
segment parameter estimate (median and 95\% credible interval) 
as a function of periodic time. Radial lines denote the start 
and end of night period if provided via the night argument.
}

\value{
  Produces a plotting frame. No value are returned.
}

\examples{
  set.seed(10)
  rate = c(rep(1,12),rep(10,12),rep(1,12),rep(10,12))
  data = rpois(48,rate)
  time = c(1:24,1:24)
  DRAW <- CircCPT(data=data, time=time)
  
  diagnostic.plot(DRAW)  
}

\keyword{diagnostic.plot}