\name{pcpt.prior}
\alias{pcpt.prior}
\title{
  Generic Function -- pcpt.prior
}
\description{
  Generic function
}
\usage{
  pcpt.prior(object)
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
  \code{\link{pcpt.prior-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
pcpt.prior(x) # retrieves the pcpt.prior slot from x
}
\keyword{ methods }
\keyword{ pcpt }
\keywoed{ internal }
