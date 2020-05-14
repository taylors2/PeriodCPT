\name{pcpt-class}
\docType{class}
\alias{pcpt-class}
\alias{pcpt-method,data.set}
\alias{pcpt-method,periodlength}
\alias{pcpt-method,minseglen}
\alias{pcpt-method,npcpts.max}
\alias{pcpt-method,distribution}
\alias{pcpt-method,pcpt.prior}
\alias{pcpt-method,param.prior}
\alias{pcpt-method,MCMC.options}
\alias{pcpt-method,MCMC.inits}
\alias{pcpt-method,MCMC.chain}
\alias{pcpt-method,MCMC.chains}
\alias{pcpt-method,n.chains}
\alias{pcpt-method,n.iter}
\alias{pcpt-method,n.burn}
\alias{pcpt-method,toggle.quiet}

\alias{pcpt-method,data.set<-}
\alias{pcpt-method,periodlength<-}
\alias{pcpt-method,minseglen<-}
\alias{pcpt-method,distribution<-}
\alias{pcpt-method,pcpt.prior<-}
\alias{pcpt-method,param.prior<-}
\alias{pcpt-method,MCMC.options<-}
\alias{pcpt-method,MCMC.inits<-}
\alias{pcpt-method,MCMC.chain<-}
\alias{pcpt-method,MCMC.chains<-}
\alias{pcpt-method,n.chains<-}
\alias{pcpt-method,n.iter<-}
\alias{pcpt-method,n.burn<-}

\title{Class "pcpt"}
\description{
  A class for within period changepoint objects.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("pcpt", ...)}.
  \describe{
    \item{\code{new("pcpt", ...)}:}{creates a new object with class pcpt }
  }
}
\section{Slots}{
  \describe{
    \item{\code{data.set}:}{Object of class \code{"ts"}, a coerced time series of the original data. }
    \item{\code{periodlength}:}{Object of class \code{"numeric"}, the period length of the data. }
    \item{\code{minseglen}:}{Object of class \code{"numeric"}, the minimum segment length (no. of observations between changepoints) used in the analysis.}
    \item{\code{npcpts.max}:}{Object of class \code{"numeric"}, maximum number of within period changepoints that can be identified. }
    \item{\code{distribution}:}{Object of class \code{"character"}, the sampling distribution of the data being analysed. }
    \item{\code{pcpt.prior}:}{Object of class \code{"list"}, containing the within period changepoint prior details.}
    \item{\code{param.prior}:}{Object of class \code{"numeric"},vector containing the within period segment parameter prior details.}
    \item{\code{MCMC.options}:}{Object of class \code{"list"}, containing algorithmic details to be supplied to the MCMC program.}
    \item{\code{MCMC.inits}:}{Object of class \code{"list"}, containing the within period changepoint positions for initiating each chain.}
    \item{\code{MCMC.chains}:}{Object of class \code{"list"}, containing the sampled MCMC chains.}
    \item{\code{param.est}:}{Object of class \code{"list"}, list where each element is a vector of parameter estimates, if requested. (TO BE IMPLEMENTED!)}
    \item{\code{pcpt.est}:}{Object of class \code{"list"}, list where each element is a vector of parameter estimates, if requested. (TO BE IMPLEMENTED!)}
    \item{\code{date}:}{Object of class \code{"character"}, date and time the changepoint analysis was run. }
        \item{\code{version}:}{Object of class \code{"character"}, version number of the package used when the analysis was run.}

  }
}

\section{Methods}{
  \describe{
    \item{data.set}{\code{signature(object = "pcpt")}: retrieves the data.set slot }
    \item{periodlength}{\code{signature(object = "pcpt")}: retrieves the periodlength slot }
    \item{minseglen}{\code{signature(object = "pcpt")}: retrieves the minseglen slot }
    \item{npcpts.max}{\code{signature(object = "pcpt")}: retrieves the npcpts.max slot }
    \item{distribution}{\code{signature(object = "pcpt")}: retrieves distribution slot }
    \item{pcpt.prior}{\code{signature(object = "pcpt")}: retrieves pcpt.prior slot }
    \item{param.prior}{\code{signature(object = "pcpt")}: retrieves param.prior slot }
    \item{MCMC.options}{\code{signature(object = "pcpt")}: retrieves MCMC.options slot }
    \item{n.chains}{\code{signature(object = "pcpt")}: retrieves from the MCMC.options slot the \code{n.chains} item. }
    \item{n.iter}{\code{signature(object = "pcpt")}: retrieves from the MCMC.options slot the \code{n.iter} item. }
    \item{n.burn}{\code{signature(object = "pcpt")}: retrieves from the MCMC.options slot the \code{n.burn} item. }
    \item{toggle.quiet}{\code{signature(object = "pcpt")}: inverts the logical item \code{quiet} in the MCMC.options slot. }
    \item{MCMC.inits}{\code{signature(object = "pcpt")}: retrieves MCMC.inits slot }
    \item{MCMC.chains}{\code{signature(object = "pcpt")}: retrieves MCMC.chains slot }
    \item{MCMC.chain}{\code{signature(object = "pcpt", index = "numeric")}: retrieves from the MCMC.chains slot the chain corresponding to the \code{index}-th chain. }

    \item{print}{\code{signature(x = "cpt")}: prints details of the pcpt object including summary}
    \item{show}{\code{signature(object = "cpt")}: prints details of the pcpt object including summary}
    \item{summary}{\code{signature(object = "pcpt", ...)}: prints a summary of the pcpt object }
    \item{plot}{\code{signature(x = "pcpt", ...)}: plots the pcpt object}

  }
}

\author{
Simon Taylor
Rebecca Killick
}

\seealso{
\code{\link{data.set-methods}},\code{\link{PeriodCPT}}
}
\examples{
showClass("pcpt") # shows the structure of the pcpt class

x = new("pcpt")       # creates a new object with the pcpt class defaults
periodlength(x)       # retrieves the periodlength slot from x
periodlength(x) <- 12 # replaces the periodlength slot from x with 12

# Example of a within period change in probability in simulated bernoulli data
set.seed(1)
x = ts(rbinom(n = 120, size = 1, rep( c(0.2,0.5), each = 6)), freq = 12)
ans = PeriodCPT(x, distribution = "bern")

print(ans) # prints details of the analysis including a summary
summary(ans)

##\/TODO\/
##########
##plot(ans) # plots the data with change (vertical line) at 100
##logLik(ans) # raw likelihood of the data with changepoints, second value is likelihood + penalty
}
\keyword{classes}
\keyword{pcpt}
\keyword{internal}
