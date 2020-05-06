\name{periodlength}
\alias{periodlength}
\title{
  Generic Function -- periodlength
}
\description{
  Generic function
}
\usage{
  periodlength(object)
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
  \code{\link{periodlength-methods}}
\examples{
x = new("pcpt") # new pcpt object
periodlength(x) # retrieves the periodlength slot from x
}
\keyword{ methods }
\keyword{ pcpt }
\keywoed{ internal }
