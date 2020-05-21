\name{summarized<-}
\alias{summarized<-}
\title{
Generic Function - summarized<-
}
\description{
Generic function
}
\usage{
summarized(object)<-value
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
\code{\link{summarized<--methods}}
}
\examples{
x=new("pcpt") # new pcpt object
summarized(x) = FALSE
summarized(x)
}

\keyword{methods}
\keyword{pcpt}
\keyword{internal}
