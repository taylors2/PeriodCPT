\name{pcpt.mode}
\alias{pcpt.mode}
\title{
  Generic Function -- pcpt.mode
}
\description{
  Generic function
}
\usage{
  pcpt.mode(object)
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
  \code{\link{pcpt.mode-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
pcpt.mode(x) = c(6, 12)
pcpt.mode(x)
}
\keyword{ methods }
\keyword{ pcpt }
\keyword{ internal }
