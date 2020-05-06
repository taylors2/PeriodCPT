\name{npcpt.max}
\alias{npcpt.max}
\title{
  Generic Function -- npcpt.max
}
\description{
  Generic function
}
\usage{
  npcpt.max(object)
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
  \code{\link{npcpt.max-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
npcpt.max(x) # retrieves the npcpt.max slot from x
}
\keyword{ methods }
\keyword{ pcpt }
\keywoed{ internal }
