param.prior.make.mean <- function(param.m, param.c, var){
  ##y ~ Norm(theta, var), theta~Norm(param.m, param.c)

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

  if(missing(var)){
    var <- 1
  }else{
    if(!is.numeric(var) || length(var) != 1 || any(var <= 0))
      stop("Known `var` parameter in distribution function is specified incorrectly.")
  }

  out <- c("var" = var, "param.m" = param.m, "param.c" = param.c)
  return(out)

}

PeriodCPT.mean <- function(data, periodlength = NULL, minseglen = 1, Mprior = c("pois", "unif"),
                           Mhyp = 1, spread = 1, param.m = 0, param.c = 1, var = 1,
                           inits = NULL, n.iter = 1e6, n.chain = 1, n.burn = 0,
                           cachesize=50, quiet=FALSE, ...){

  distribution <- "mean"
  if(!is.numeric(data))
    stop("Data is invalid for Normal sampling distribution.")

  Mprior <- match.arg(Mprior)
  ans <- class_input(data = data, periodlength = periodlength, minseglen = minseglen,
                     distribution = distribution, Mprior = Mprior, Mhyp = Mhyp,
                     spread = spread, inits = inits,
                     n.iter, n.chain, n.burn, cachesize, quiet,
                     param.m, param.c, var, ...)

  ans <- PeriodCPT.main(ans)
  ans <- eval(paste0("SummariseOutput.",distribution,"(ans)"))
  return(ans)
}


SummariseOutput.mean <- function(object){
  return(object)
}

