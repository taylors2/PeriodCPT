\name{seglen}
\alias{seglen}
\title{
  Generic Function -- seglen
}
\description{
  Generic function
}
\usage{
  seglen(object)
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
  \code{\link{seglen-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
periodlength(x) = 12
pcpt.mode(x) = c(6, 12)
seglen(x)
}
\keyword{ methods }
\keyword{ pcpt }
\keyword{ internal }
