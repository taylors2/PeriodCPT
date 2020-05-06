\name{param.prior}
\alias{param.prior}
\title{
  Generic Function -- param.prior
}
\description{
  Generic function
}
\usage{
  param.prior(object)
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
  \code{\link{param.prior-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
param.prior(x) # retrieves the param.prior slot from x
}
\keyword{ methods }
\keyword{ pcpt }
\keywoed{ internal }
