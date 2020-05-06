\name{n.burn}
\alias{n.burn}
\title{
  Generic Function -- n.burn
}
\description{
  Generic function
}
\usage{
  n.burn(object)
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
  \code{\link{n.burn-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
MCMC.options(x) = MCMC.options.make()
n.burn(x) # retrieves the n.burn item in the MCMC.options slot from x
}
\keyword{ methods }
\keyword{ pcpt }
\keywoed{ internal }
