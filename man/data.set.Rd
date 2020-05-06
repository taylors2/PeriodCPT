\name{data.set}
\alias{data.set}
\title{
  Generic Function -- data.set
}
\description{
  Generic function
}
\usage{
  data.set(object)
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
  \code{\link{data.set-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
data.set(x) # retrieves the data.set slot from x
}
\keyword{ methods }
\keyword{ pcpt }
\keywoed{ internal }
