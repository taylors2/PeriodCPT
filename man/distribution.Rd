\name{distribution}
\alias{distribution}
\title{
  Generic Function -- distribution
}
\description{
  Generic function
}
\usage{
  distribution(object)
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
  \code{\link{distribution-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
distribution(x) # retrieves the distribution slot from x
}
\keyword{ methods }
\keyword{ pcpt }
\keywoed{ internal }
