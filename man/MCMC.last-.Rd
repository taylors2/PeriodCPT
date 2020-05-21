\name{MCMC.last<-}
\alias{MCMC.last<-}
\title{
Generic Function - MCMC.last<-
}
\description{
Generic function
}
\usage{
MCMC.last(object)<-value
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
\code{\link{MCMC.last<--methods}}
}
\examples{
x=new("pcpt") # new pcpt object
MCMC.last(x) = list(`1` = 1, `2` = c(6, 12))
MCMC.last(x)
}

\keyword{methods}
\keyword{pcpt}
\keyword{internal}
