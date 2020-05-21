\name{summarised}
\alias{summarised}
\title{
  Generic Function -- summarised
}
\description{
  Generic function
}
\usage{
  summarised(object)
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
  \code{\link{summarised-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
summarised(x) = FALSE
summarised(x)
}
\keyword{ methods }
\keyword{ pcpt }
\keyword{ internal }
