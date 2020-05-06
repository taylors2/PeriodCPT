param.prior.make.norm <- function(param.m, param.c, param.a, param.b, ...){
  ##y ~ Norm(theta[1], theta[2]),
  ##theta[1] ~ Norm(param.m, param.c * theta[2])
  ##theta[2] ~ IGamma(param.a, param.b)

  if(missing(param.m)){
    param.m <- 0
  }else{
    if(!is.numeric(param.m) || length(param.m) != 1)
      stop("Normal-InvGamma `param.m` hyper parameter specified incorrectly.")
  }

  if(missing(param.c)){
    param.c <- 1
  }else{
    if(!is.numeric(param.c) || length(param.c) != 1 || any(param.c <= 0))
      stop("Normal-InvGamma `param.c` hyper parameter specified incorrectly.")
  }

  if(missing(param.a)){
    param.a <- 1
  }else{
    if(!is.numeric(param.a) || length(param.a) != 1 || any(param.a <= 0))
      stop("Normal-InvGamma `param.a` hyper parameter specified incorrectly.")
  }

  if(missing(param.b)){
    param.b <- 1
  }else{
    if(!is.numeric(param.b) || length(param.b) != 1 || any(param.b <= 0))
      stop("Normal-InvGamma `param.b` hyper parameter specified incorrectly.")
  }

  out <- c("param.m" = param.m, "param.c" = param.c,
           "param.a" = param.a, "param.b" = param.b)
  return(out)
}


PeriodCPT.norm <- function(data, periodlength = NULL, minseglen = 1, Mprior = c("pois", "unif"),
                           Mhyp = 1, spread = 1, param.m = 0, param.c = 1, param.a = 1,
                           param.b = 1, inits = NULL, n.iter = 1e6, n.chain = 1, n.burn = 0,
                           cachesize=50, quiet=FALSE, ...){

  distribution <- "norm"
  if(!is.numeric(data))
    stop("Data is invalid for Normal sampling distribution.")

  Mprior <- match.arg(Mprior)
  ans <- class_input(data = data, periodlength = periodlength, minseglen = minseglen,
                     distribution = distribution, Mprior = Mprior, Mhyp = Mhyp,
                     spread = spread, inits = inits,
                     n.iter, n.chain, n.burn, cachesize, quiet,
                     param.m, param.c, param.a, param.b, ...)

  ans <- PeriodCPT.main(ans)
  ans <- eval(paste0("SummariseOutput.",distribution,"(ans)"))
  return(ans)
}

SummariseOutput.norm <- function(object){
  return(object)
}

