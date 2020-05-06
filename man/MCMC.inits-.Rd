\name{MCMC.inits<-}
\alias{MCMC.inits<-}
\title{
Generic Function - MCMC.inits<-
}
\description{
Generic function
}
\usage{
MCMC.inits(object)<-value
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
\code{\link{MCMC.inits<--methods}}
}
\examples{
x=new("pcpt") # new pcpt object
MCMC.inits(x) = Definie.inits(x) # replaces the existing MCMC.inits slot in x
}

\keyword{methods}
\keyword{pcpt}
\keyword{internal}
