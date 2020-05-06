\name{MCMC.chain}
\alias{MCMC.chain}
\title{
  Generic Function -- MCMC.chain
}
\description{
  Generic function
}
\usage{
  MCMC.chain(object, index)
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
  \code{\link{MCMC.chain-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
MCMC.chains(x) = list(`1` = NULL, `2` = NULL)
MCMC.chain(x, 2) # retrieves the the second listed item from the MCMC.chains slot from x
}
\keyword{ methods }
\keyword{ pcpt }
\keyword{ internal }
