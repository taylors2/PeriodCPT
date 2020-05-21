\name{param.mode}
\alias{param.mode}
\title{
  Generic Function -- param.mode
}
\description{
  Generic function
}
\usage{
  param.mode(object)
}
\arguments{
  \item{object}{
    Depending on the class of \code{object} depends on the method used (and if one exists)
  }
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
  \code{\link{param.mode-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
param.mode(x) = matrix(c(0.2, 0.8), nrow = 1, ncol = 2)
param.mode(x)
}
\keyword{ methods }
\keyword{ pcpt }
\keyword{ internal }
