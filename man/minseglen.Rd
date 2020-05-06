\name{minseglen}
\alias{minseglen}
\title{
  Generic Function -- minseglen
}
\description{
  Generic function
}
\usage{
  minseglen(object)
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
  \code{\link{minseglen-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
minseglen(x) # retrieves the minseglen slot from x
}
\keyword{ methods }
\keyword{ pcpt }
\keyword{ internal }
