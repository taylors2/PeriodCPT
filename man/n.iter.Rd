\name{n.iter}
\alias{n.iter}
\title{
  Generic Function -- n.iter
}
\description{
  Generic function
}
\usage{
  n.iter(object)
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
  \code{\link{n.iter-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
MCMC.options(x) = MCMC.options.make()
n.iter(x) # retrieves the n.iter item in the MCMC.options slot from x
}
\keyword{ methods }
\keyword{ pcpt }
\keyword{ internal }
