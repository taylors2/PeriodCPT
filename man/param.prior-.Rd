\name{param.prior<-}
\alias{param.prior<-}
\title{
Generic Function - param.prior<-
}
\description{
Generic function
}
\usage{
param.prior(object)<-value
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
\code{\link{param.prior<--methods}}
}
\examples{
x=new("pcpt") # new pcpt object
param.prior(x) = c("param.a" = 1, "param.b" = 1)
   # replaces the existing param.prior slot in x
}

\keyword{methods}
\keyword{pcpt}
\keyword{internal}
