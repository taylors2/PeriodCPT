\name{MCMC.options<-}
\alias{MCMC.options<-}
\title{
Generic Function - MCMC.options<-
}
\description{
Generic function
}
\usage{
MCMC.options(object)<-value
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
\code{\link{MCMC.options<--methods}}
}
\examples{
x=new("pcpt") # new pcpt object
MCMC.options(x) = MCMC.options.make()
   # replaces the existing MCMC.options slot in x
}

\keyword{methods}
\keyword{pcpt}
\keyword{internal}
