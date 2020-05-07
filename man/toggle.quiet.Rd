\name{toggle.quiet}
\alias{toggle.quiet}
\title{
  Generic Function -- toggle.quiet
}
\description{
  Generic function
}
\usage{
  toggle.quiet(object)
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
  \code{\link{toggle.quiet-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
MCMC.options(x) = list(
         n.iter = 1e6, n.chain = 1,
         n.burn = 1e4, cachesize = 50,
         quiet = FALSE)
x@MCMC.options$quiet
toggle.quiet(x) # inverts the logical item quiet in the MCMC.options slot.
x@MCMC.options$quiet
}
\keyword{ methods }
\keyword{ pcpt }
\keyword{ internal }
