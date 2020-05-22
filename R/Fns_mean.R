param.prior.make.mean <- function(param.m, param.c, const.var, ...){
  ##y ~ Norm(theta, const.var), theta~Norm(param.m, param.c)

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

  if(missing(const.var)){
    const.var <- 1
  }else{
    if(!is.numeric(const.var) || length(const.var) != 1 || any(const.var <= 0))
      stop("Known `const.var` parameter in distribution function is specified incorrectly.")
  }

  out <- c("param.m" = param.m, "param.c" = param.c, "const.var" = const.var)
  return(out)

}

PeriodCPT.mean <- function(data, periodlength = NULL, minseglen = 1, Mprior = c("pois", "unif"),
                           Mhyp = 1, spread = 1, param.m = 0, param.c = 1, const.var = 1,
                           inits = NULL, n.iter = 1e6, n.chains = 1, n.burn = 0,
                           cachesize=50, quiet=FALSE, ...){

  distribution <- "mean"
  nsegparam <- 1
  if(!is.numeric(data))
    stop("Data is invalid for Normal sampling distribution.")

  Mprior <- match.arg(Mprior)

  ans <- class_input(data = data, periodlength = periodlength, minseglen = minseglen,
                     distribution = distribution, nsegparam = nsegparam, Mprior = Mprior, Mhyp = Mhyp,
                     spread = spread, inits = inits, n.iter = n.iter,
                     n.chains = n.chains, n.burn = n.burn, cachesize = cachesize,
                     quiet = quiet, param.m = param.m, param.c = param.c,
                     const.var = const.var, ...)

  ans <- PeriodCPT.main(ans)
  return(ans)
}

