\name{result}
\alias{result}
\title{
  Generic Function -- result
}
\description{
  Generic function
}
\usage{
  result(object, index)
}
\arguments{
  \item{object}{
    Depending on the class of \code{object} depends on the method used (and if one exists)
  }
  \item{index}{ Single positive integer giving chain index.}
}
\details{
  Generic Function
}
\value{
  Depends on the calss of \code{object}, see individual methods
}
\author{
  Simon Taylor
  Rebecca Killick
}
\seealso{
  \code{\link{result-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
results(x) = list(`1` = NULL, `2` = NULL)
result(x, 2) # retrieves the the second listed item from the results slot from x
}
\keyword{ methods }
\keyword{ pcpt }
\keyword{ internal }
