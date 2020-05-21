\name{summarized}
\alias{summarized}
\title{
  Generic Function -- summarized
}
\description{
  Generic function
}
\usage{
  summarized(object)
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
  \code{\link{summarized-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
summarized(x) = FALSE
summarized(x)
}
\keyword{ methods }
\keyword{ pcpt }
\keyword{ internal }
