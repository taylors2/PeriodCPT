\name{pcpt.prior<-}
\alias{pcpt.prior<-}
\title{
Generic Function - pcpt.prior<-
}
\description{
Generic function
}
\usage{
pcpt.prior(object)<-value
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
\code{\link{pcpt.prior<--methods}}
}
\examples{
x=new("pcpt") # new pcpt object
pcpt.prior(x) = list(Mprior = "pois", Mhyp = 1, spread = 1)
   # replaces the existing pcpt.prior slot in x
}

\keyword{methods}
\keyword{pcpt}
\keyword{internal}
