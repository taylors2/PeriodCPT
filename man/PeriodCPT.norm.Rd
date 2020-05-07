\name{PeriodCPT.norm}
\alias{PeriodCPT.norm}
\title{
  ~~function to evaluate within period changepoint events for normal data~~
}
\description{
  Evaluation of within period changepoint events for Normal distributed observations where interest lies in both mean and variance.
}
\usage{
PeriodCPT.norm(data, periodlength = NULL, minseglen = 1, Mprior = c("pois", "unif"),
                           Mhyp = 1, spread = 1, param.m = 0, param.c = 1, param.a = 1,
                           param.b = 1, inits = NULL, n.iter = 1e6, n.chains = 1,
                           n.burn = 0, cachesize = 50, quiet = FALSE, ...)
}

\arguments{
  \item{data}{
    A vector or ts object containing the data within which you which to find within period changepoint events.
  }
  \item{periodlength}{
    Positive integer specifying the period length. If NULL, then information is obtained from ts object specified by data.
  }
  \item{minseglen}{Positive integer giving the minumum segment length (no. of observations between within-period changes), set at 1 for default.}
  \item{Mprior}{Option specifying the prior for the number of within period changepoint events.}
  \item{Mhyp}{Numeric vector giving the hyper-parameters for the number of within period chagepoint prior, see details.}
  \item{spread}{Positive number specifying the prior degree of spread around the period for the changepoint positions. By default, 1 indicates random uniform spread.}
  \item{param.m}{A single number giving the location hyper-parameters for the normal-InvGamma prior distribution on the within period segment parameter. Set to 0 by default.}
  \item{param.c, param.a, param.b}{Positive numbers giving the scale & shape hyper-parameters for the normal-InvGamma prior distribution on the within period segment parameter. All set to 1 by default.}
  \item{inits}{Initial values for the chains. See details}
  \item{n.iter}{Positive integer giving the number of iterations for the MCMC.}
  \item{n.chains}{Positive integer giving the number of chains to evaluate.}
  \item{n.burn}{Positive integer giving the initial number iterations to discard as the burn-in period.}
  \item{cachesize}{Positive integer giving the cache size to minimise re-evaluation of repeated scenarios within the MCMC algorithm. Set to 50 by default.}
  \item{quiet}{Logical flag indicating whether a progress bar should be printed}
  \item{...}{Additional inputs, see details.}
}

\details{
  This function is used to find with period changepoint events for selected sampling distribution and statistic. A changepoint is denoted as the (!)LAST OBSERVATION OF THE CURRENT SEGMENT/REGIME(!)

  Segment parameters: This model assumes that, for a given segment, data~norm(mu, var) with prior distribution (mu, var)~normal-InvGamma(prior.m, prior.c, prior.a, prior.b).

  Within period changepoint parameters: The prior for the number of within period changepoint events is either:
  \describe{
  \item{pois}{
    Poisson distribution truncated to be between 1 and the maximum number permitted within the period given the minimum segment length constraint. In this case, \code{Mhyp} should contain a single positive number that denotes the lambda hyper-parameter.
  }
  \item{unif}{
    Discrete uniform distribution on the integers from 1 to the maximum number permitted within the period given the minimum segment length constraint. The argument \code{Mhyp} is ignored.
  }
  }
  Note: the scenario where the number of within period changepoints is equal to 1 defines the case where there are no within-period changepoint events and so model reduces to just the sampling distribution.
  Given the total number of changepoint events, the prior for their positions is given in terms of their excess segment length beyond the minimum condition: delta[i] = pcpt[i] - pcpt[i-1] - minseglen. The prior on the vector of excess segment lengths is a Dirichlet-multinomial distribution with shape vector given by \code{spread} (same size for all components). If \code{spread}<1 then there is a prior preference that the changepoints cluster, \code{spread}>1 indicates a prior preference for repellent events, else \code{spread}=1 specifies a uniform distribution. The excess segment lengths are mapped back to the within period changepoint positions with the aid of a anchoring event that follows a uniform distribution over the within peroid positions. If necessary, events are permuted such that all are defined within the interval [0, periodlength).

  The chain inital values, \code{init}, can be specified in one of the following ways:
  \describe{
  \item{NULL}{Default, positions are generated from the prior for each chain.}
  \item{"ends"}{This option specifies inital values for two chains (adjusts \code{n.chains()} accordingly) where the first denotes the trivial case of no within period changepoint events, and the other denots a packed case by fitting in as many events as possible with respect to the minimum segment length.}
  \item{vector}{A vector of integers giving the inital value for one (and only one) requested chain.}
  \item{list}{A list where each item specifies the vector of initial values for each chain.}
  \item{function}{The initalisation of a given chain is determined/generated by a provided function. Input arguments for this function can be specified via \code{...} of the main PeriodCPT funciton. The input names \code{pcpt.object} and \code{chain.index} can be specified to respectively obtain information held within the internally generated \code{pcpt} class object or to specify a function that depends on a specific chain via an integer index.
  }
  }
  Note: For codeing ease, the trivial case of no changepoint events is specified by \code{c(0)}.
}

\value{
  An object of S4 class "pcpt" is returned.
}

\author{
Simon Taylor
}

\seealso{
  \code{\linkS4class{pcpt}}, \code{\link{PeriodCPT}}
}

\examples{
  data = ts( rnorm(120, mean = rep(c(-5, 5), each=6),
    sd = rep(c(1, 3), each=6)), freq = 1/12)
  ans  = PeriodCPT.norm(data)
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
\keyword{changepoint}
