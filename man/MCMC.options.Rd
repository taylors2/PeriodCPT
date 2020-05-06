\name{MCMC.options}
\alias{MCMC.options}
\title{
  Generic Function -- MCMC.options
}
\description{
  Generic function
}
\usage{
  MCMC.options(object)
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
  \code{\link{MCMC.options-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
MCMC.options(x) # retrieves the MCMC.options slot from x
}
\keyword{ methods }
\keyword{ pcpt }
\keyword{ internal }
