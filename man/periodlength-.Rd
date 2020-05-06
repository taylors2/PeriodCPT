\name{periodlength<-}
\alias{periodlength<-}
\title{
Generic Function - periodlength<-
}
\description{
Generic function
}
\usage{
periodlength(object)<-value
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
\code{\link{periodlength<--methods}}
}
\examples{
x=new("pcpt") # new pcpt object
periodlength(x) <- 12 # replaces the existing periodlength slot in x with 12
}

\keyword{methods}
\keyword{pcpt}
\keyword{internal}
