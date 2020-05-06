\name{MCMC.chain<-}
\alias{MCMC.chain<-}
\title{
Generic Function - MCMC.chain<-
}
\description{
Generic function
}
\usage{
MCMC.chain(object)<-value
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
\code{\link{MCMC.chain<--methods}}
}
\examples{
x=new("pcpt") # new pcpt object
MCMC.chains(x) = list(`1` = NULL, `2` = NULL)
MCMC.chain(x, 2) = 0  #assigns the number zero to the second list item in MCMC.chains slot of x
}

\keyword{methods}
\keyword{pcpt}
\keyword{internal}
