\name{quantile-methods}
\docType{methods}
\alias{quantile-methods}
\alias{quantile,pcpt-method}
\title{ ~~ Methods for Function print in Package `base' ~~}
\description{
 ~~ Methods for function \code{quantile} in Package `base' ~~
}
\section{Methods}{
\describe{

\item{\code{signature(x = "pcpt", ...)}}{
	Evaluates the requires quantile estimate for the segment parameters for each within period time step.
}
}}

\value{
Returns a matrix with each row containing the estimated quantiles (according to the requested probabilities specified as a numerical vector using the optional \code{prob} argument) of the marginal posterior distribution for each of the segment parameters for each time step within the period length interval.
}

\examples{
set.seed(1)
x  = ts(rbinom(240, size = 1, prob = rep(rep( c(0.2, 0.8), each = 12), times = 10)), frequency = 24)
ans = PeriodCPT(data = x, distribution = "bern", quiet = TRUE, n.iter = 1000)
quantile(ans)
}

\keyword{methods}
\keyword{quantile}
\keyword{changepoint}
\keyword{internal}
