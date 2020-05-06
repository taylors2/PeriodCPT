\name{MCMC.chains}
\alias{MCMC.chains}
\title{
  Generic Function -- MCMC.chains
}
\description{
  Generic function
}
\usage{
  MCMC.chains(object)
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
  \code{\link{MCMC.chains-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
MCMC.chains(x) # retrieves the MCMC.chains slot from x
}
\keyword{ methods }
\keyword{ pcpt }
\keywoed{ internal }
