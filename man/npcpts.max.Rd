\name{npcpts.max}
\alias{npcpts.max}
\title{
  Generic Function -- npcpts.max
}
\description{
  Generic function
}
\usage{
  npcpts.max(object)
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
  \code{\link{npcpts.max-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
npcpts.max(x) # retrieves the npcpts.max slot from x
}
\keyword{ methods }
\keyword{ pcpt }
\keyword{ internal }
