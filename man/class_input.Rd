\name{class_input}
\alias{class_input}
\title{
Input all required arguments into pcpt class - Only intended for developer use.
}
\description{
This function helps to input all the necessary information into the correct format for the \code{pcpt} class.

This function is called by \code{PeriodCPT} and its distribution specific functions. This is not intended for use by regular users of the package.  It is exported for developers to call directly for speed and convenience.

WARNING: Checks on inputs may not be comprehensive!
}
\usage{
class_input <- function(data, periodlength, minseglen, distribution,
                        Mprior, Mhyp, spread, inits, ...)
}
\arguments{

  \item{data}{
    A vector or ts object containing the data within which you which to find within period changepoint events. If not specified, the period length is taken from the ts object.
  }
  \item{periodlength}{
    Positive integer specifying the period length. Information is obtained from data ts-object if not specified.
  }

  \item{minseglen}{Positive integer giving the minumum segment length (no. of observations between within-period changes).}

  \item{distribution}{
    Character specifying the sampling distribution/statistic to investigate within period changepoint events.
  }

  \item{Mprior}{Character specifying the prior for the number of within period changepoint events.}
  \item{Mhyp}{Numeric vector giving the hyper-parameters for the number of within period chagepoint prior.}
  \item{spread}{Positive number specifying the prior degree of spread around the period for the changepoint positions.}
  \item{inits}{Initial values for the chains.}
  \item{...}{
    Additional optional arguments that specifies prior details, MCMC parameters and output format.
  }


}


\details{
	This function takes all the input required for the \code{pcpt} class and enters it into the object.

		This function is exported for developer use only.  It perform some checks on inputs, but it is simply a convenience function for converting the output of the worker functions into a nice format for the \code{pcpt} class.
}
\value{
  An object of class \code{pcpt} filled with the given attributes.
}
\author{
Simon Taylor
Rebecca Killick
}
\seealso{
\code{\link{PeriodCPT}},\code{\linkS4class{pcpt}}
}
\examples{
#This function should only be used by developers, see its use in PeriodCPT and distribution specific variations.
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
\keyword{changepoint}
