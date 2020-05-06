\name{minseglen<-}
\alias{minseglen<-}
\title{
Generic Function - minseglen<-
}
\description{
Generic function
}
\usage{
minseglen(object)<-value
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
\code{\link{minseglen<--methods}}
}
\examples{
x=new("pcpt") # new pcpt object
minseglen(x) <- 4 # replaces the existing minseglen slot in x with 4
}

\keyword{methods}
\keyword{pcpt}
\keyword{internal}
