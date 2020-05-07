param.prior.make.bern <- function(param.a, param.b, ...){
  ##y~bern(theta), theta~beta(param.a, param.b)
  if(missing(param.a)){
    param.a <- 1
  }else{
    if(!is.numeric(param.a) || length(param.a) != 1 || any(param.a <= 0))
      stop("Beta `param.a` hyper parameter specified incorrectly.")
  }

  if(missing(param.b)){
    param.b <- 1
  }else{
    if(!is.numeric(param.b) || length(param.b) != 1 || any(param.b <= 0))
      stop("Beta `param.b` hyper parameter specified incorrectly.")
  }
  out <- c("param.a" = param.a, "param.b" = param.b)
  return(out)
}


PeriodCPT.bern <- function(data, periodlength = NULL, minseglen = 1, Mprior = c("pois", "unif"),
                           Mhyp = 1, spread = 1, param.a = 1, param.b = 1,
                           inits = NULL, n.iter = 1e6, n.chains = 1, n.burn = 0,
                           cachesize = 50, quiet = FALSE, ...){

  distribution <- "bern"
  if(!is.numeric(data))
    stop("Data is invalid for Bernoulli sampling distribution.")
  if(!all(data == 0L | data == 1L))
    stop("Data is invalid for Bernoulli sampling distribution.")

  Mprior <- match.arg(Mprior)
  ans <- class_input(data = data, periodlength = periodlength, minseglen = minseglen,
                     distribution = distribution, Mprior = Mprior, Mhyp = Mhyp,
                     spread = spread, inits = inits, n.iter = n.iter, n.chains = n.chains,
                     n.burn = n.burn, cachesize = cachesize, quiet = quiet,
                     param.a = param.a, param.b = param.b, ...)

  ans <- PeriodCPT.main(ans)
  return(ans)
}

SummariseOutput.bern <- function(object){
  return(object)
}

