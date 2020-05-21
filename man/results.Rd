\name{results}
\alias{results}
\title{
  Generic Function -- results
}
\description{
  Generic function
}
\usage{
  results(object)
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
  \code{\link{results-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
results(x) = list(`1` = NULL, `2` = NULL)
result(x, 2) # retrieves the the second listed item from the MCMC.chains slot from x
}
\keyword{ methods }
\keyword{ pcpt }
\keyword{ internal }
