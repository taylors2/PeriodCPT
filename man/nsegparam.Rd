\name{nsegparam}
\alias{nsegparam}
\title{
  Generic Function -- nsegparam
}
\description{
  Generic function
}
\usage{
  nsegparam(object)
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
  \code{\link{nsegparam-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
nsegparam(x) = 1
nsegparam(x)
}
\keyword{ methods }
\keyword{ pcpt }
\keyword{ internal }
