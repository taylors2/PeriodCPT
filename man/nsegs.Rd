\name{nsegs}
\alias{nsegs}
\title{
  Generic Function -- nsegs
}
\description{
  Generic function
}
\usage{
  nsegs(object)
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
  \code{\link{nsegs-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
periodlength(x) = 12
pcpt.mode(x) = c(6, 12)
nsegs(x)
}
\keyword{ methods }
\keyword{ pcpt }
\keyword{ internal }
