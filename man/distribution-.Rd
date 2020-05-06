\name{distribution<-}
\alias{distribution<-}
\title{
Generic Function - distribution<-
}
\description{
Generic function
}
\usage{
distribution(object)<-value
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
\code{\link{distribution<--methods}}
}
\examples{
x=new("pcpt") # new pcpt object
distribution(x) <- "bern" # replaces the existing distribution slot in x with "bern"
}

\keyword{methods}
\keyword{pcpt}
\keyword{internal}
