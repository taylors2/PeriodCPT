\name{quiet}
\alias{quiet}
\title{
  Generic Function -- quiet
}
\description{
  Generic function
}
\usage{
  quiet(object)
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
  \code{\link{quiet-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
MCMC.options(x) = list(
         n.iter = 1e6, n.chains = 1,
         n.burn = 1e4, cachesize = 50,
         quiet = FALSE)
quiet(x) # retrieves the quiet item in the MCMC.options slot from x
}
\keyword{ methods }
\keyword{ pcpt }
\keyword{ internal }
