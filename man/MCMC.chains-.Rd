\name{MCMC.chains<-}
\alias{MCMC.chains<-}
\title{
Generic Function - MCMC.chains<-
}
\description{
Generic function
}
\usage{
MCMC.chains(object)<-value
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
\code{\link{MCMC.chains<--methods}}
}
\examples{
x=new("pcpt") # new pcpt object
MCMC.chains(x) = vector(`1` = NULL, `2` = NULL)
    # replaces the existing MCMC.chains slot in x
}

\keyword{methods}
\keyword{pcpt}
\keyword{internal}
