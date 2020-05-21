\name{fit<-}
\alias{fit<-}
\title{
Generic Function - fit<-
}
\description{
Generic function
}
\usage{
fit(object)<-value
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
\code{\link{fit<--methods}}
}
\examples{
x=new("pcpt") # new pcpt object
fit(x) = 100
fit(x)
}

\keyword{methods}
\keyword{pcpt}
\keyword{internal}
