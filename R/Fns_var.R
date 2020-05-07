param.prior.make.var <- function(param.a, param.b, ...){
  ##y ~ Norm(0, theta), theta ~ InvGama(param.a, param.b)

  if(missing(param.a)){
    param.a <- 1
  }else{
    if(!is.numeric(param.a) || length(param.a) != 1 || any(param.a <= 0))
      stop("InvGamma `param.a` hyper parameter specified incorrectly.")
  }

  if(missing(param.b)){
    param.b <- 1
  }else{
    if(!is.numeric(param.b) || length(param.b) != 1 || any(param.b <= 0))
      stop("InvGamma `param.b` hyper parameter specified incorrectly.")
  }
  out <- c("param.a" = param.a, "param.b" = param.b)
  return(out)
}


PeriodCPT.var <- function(data, periodlength = NULL, minseglen = 1, Mprior = c("pois", "unif"),
                          Mhyp = 1, spread = 1, param.a = 1, param.b = 1,
                          inits = NULL, n.iter = 1e6, n.chains = 1, n.burn = 0,
                          cachesize=50, quiet=FALSE, ...){

  distribution <- "var"
  if(!is.numeric(data))
    stop("Data is invalid for Normal sampling distribution.")

  Mprior <- match.arg(Mprior)
  ans <- class_input(data = data, periodlength = periodlength, minseglen = minseglen,
                     distribution = distribution, Mprior = Mprior, Mhyp = Mhyp,
                     spread = spread, inits = inits, n.iter = n.iter, n.chains = n.chains,
                     n.burn = n.burn, cachesize = cachesize, quiet = quiet,
                     param.a = param.a, param.b = param.b, ...)

  ans <- PeriodCPT.main(ans)
  return(ans)
}


SummariseOutput.var <- function(object){
  return(object)
}
