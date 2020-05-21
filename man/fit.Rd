\name{fit}
\alias{fit}
\title{
  Generic Function -- fit
}
\description{
  Generic function
}
\usage{
  fit(object)
}
\arguments{
  \item{object}{
    Depending on the class of \code{object} depends on the method used (and if one exists)
  }
}
\details{
  Generic Function

  The five fit measurements defined below are based on the most common sampled within period changepoint vector:
  \describe{
  \item{1}{ Log of the marginalised sampling distribution and segment prior distributions. This is used in the core of the PeriodCPT algorithm.}
  \item{2}{ Log of maximum value from the sampling distribution and segment prior distributions. This point corresponds to the posterior mode of the segment parameters.}
  \item{3}{ Log of the integrated sampling distribution over the segment parameters.}
  \item{4}{ Log of the sampling distribution evaluated at the segment parameter posterior mode.}
  \item{5}{ Log of the maximum value of the sampling distribution. This corresponds to evaluation at the maximum likelihood values.}
  }
  All the above measures are multiplied by -2 to give the returned fit values.


}
\value{
  Depends on the calss of \code{object}, see individual methods
}
\author{
  Simon Taylor
  Rebecca Killick
}
\seealso{
  \code{\link{fit-methods}}
}
\examples{
x = new("pcpt") # new pcpt object
fit(x) = 100
fit(x)
}
\keyword{ methods }
\keyword{ pcpt }
\keyword{ internal }
