\name{n.iter<-}
\alias{n.iter<-}
\title{
Generic Function - n.iter<-
}
\description{
Generic function
}
\usage{
n.iter(object)<-value
}
\arguments{
  \item{object}{
	Depending on the class of \code{object} depends on the method used (and if one exists)
}
  \item{value}{
	Replacement value
}
}
\details{
	Generic Function
}
\value{
	Depends on the class of \code{object}, see individual methods
}
\author{
Simon Taylor
Rebecca Killick
}


\seealso{
\code{\link{n.iter<--methods}}
}
\examples{
x=new("pcpt") # new pcpt object
MCMC.options(x) = list(
         n.iter = 1e6, n.chains = 1,
         n.burn = 1e4, cachesize = 50,
         quiet = FALSE)
n.iter(x) = 1e7 # replaces the existing n.iter item in the MCMC.options slot in x
}

\keyword{methods}
\keyword{pcpt}
\keyword{internal}
