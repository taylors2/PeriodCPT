\name{result<-}
\alias{result<-}
\title{
Generic Function - result<-
}
\description{
Generic function
}
\usage{
result(object, index)<-value
}
\arguments{
  \item{object}{
	Depending on the class of \code{object} depends on the method used (and if one exists)
	}
	\item{index}{ Single positive integer giving chain index.}

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
\code{\link{result<--methods}}
}
\examples{
x=new("pcpt") # new pcpt object
results(x) = list(`1` = NULL, `2` = NULL)
result(x, 2) = 1  #assigns the number 1 to the second list item in results slot of x
}

\keyword{methods}
\keyword{pcpt}
\keyword{internal}
