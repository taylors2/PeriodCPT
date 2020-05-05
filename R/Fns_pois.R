param.prior.make.pois <- function(){
  ##y~pois(lambda), lambda~gamma(param.a, param.b)
  if(missing(param.a)){
    param.a <- 1
  }else{
    if(!is.numeric(param.a) || length(param.a) != 1 || any(param.a <= 0))
      stop("Gamma `param.a` hyper parameter specified incorrectly.")
  }

  if(missing(param.b)){
    param.b <- 1
  }else{
    if(!is.numeric(param.b) || length(param.b) != 1 || any(param.b <= 0))
      stop("Gamma `param.b` hyper parameter specified incorrectly.")
  }
  out <- c("param.a" = param.a, "param.b" = param.b)
  return(out)
}

PeriodCPT.pois <- function(data, periodlength = NULL, minseglen = 1, Mprior = c("pois", "unif"),
                           Mhyp = 1, spread = 1, param.a = 1, param.b = 1,
                           inits = NULL, n.iter = 1e6, n.chain = 1, n.burn = 0,
                           cachesize=50, quiet=FALSE){

  distribution <- "pois"
  if(!is.numeric(data))
    stop("Data is invalid for Poisson sampling distribution.")
  if(any(data < 0))
    stop("Data is invalid for Poisson sampling distribution.")
  if(any(data != floor(data)))
    stop("Data is invalid for Poisson sampling distribution.")

  Mprior <- match.arg(Mprior)
  ans <- class_input(data = data, periodlength = periodlength, minseglen = minseglen,
                     distribution = distribution, Mprior = Mprior, Mhyp = Mhyp,
                     spread = spread, inits = inits,
                     n.iter, n.chain, n.burn, cashesize, quiet,
                     param.a, param.b)

  ans <- PeriodCPT.main(ans)
  ans <- eval(paste0("SummariseOutput.",distribution,"(ans)"))
  return(ans)
}

SummariseOutput.pois <- function(object){
  return(object)
}
