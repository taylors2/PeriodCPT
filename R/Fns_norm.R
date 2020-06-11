param.prior.make.norm <- function(param.m, param.c, param.a, param.b, ...){
  ##y ~ Norm(theta[1], theta[2]),
  ##theta[1] ~ Norm(param.m, param.c * theta[2])
  ##theta[2] ~ IGamma(param.a, param.b)

  if(missing(param.m)){
    param.m <- 0
  }else{
    if(!is.numeric(param.m) || length(param.m) != 1)
      stop("Hyper-parameter `param.m` specified incorrectly.")
  }

  if(missing(param.c)){
    param.c <- 1
  }else{
    if(!is.numeric(param.c) || length(param.c) != 1 || any(param.c <= 0))
      stop("Hyper-parameter `param.c` specified incorrectly.")
  }

  if(missing(param.a)){
    param.a <- 1
  }else{
    if(!is.numeric(param.a) || length(param.a) != 1 || any(param.a <= 0))
      stop("Hyper-parameter `param.a` specified incorrectly.")
  }

  if(missing(param.b)){
    param.b <- 1
  }else{
    if(!is.numeric(param.b) || length(param.b) != 1 || any(param.b <= 0))
      stop("Hyper-parameter `param.b` specified incorrectly.")
  }

  out <- c("param.m" = param.m, "param.c" = param.c,
           "param.a" = param.a, "param.b" = param.b)
  return(out)
}


PeriodCPT.norm <- function(data, periodlength = NULL, minseglen = 1, Mprior = c("pois", "unif"),
                           Mhyp = 1, spread = 1, param.m = 0, param.c = 1, param.a = 1,
                           param.b = 1, inits = NULL, n.iter = 1e6, n.chains = 1, n.burn = 0,
                           cachesize=50, quiet=FALSE, ...){

  distribution <- "norm"
  nsegparam <- 2
  Mprior <- Mprior[1]
  ans <- class_input(data = data, periodlength = periodlength, minseglen = minseglen,
                     distribution = distribution, nsegparam = nsegparam, Mprior = Mprior, Mhyp = Mhyp,
                     spread = spread, inits = inits, n.iter = n.iter, n.chains = n.chains,
                     n.burn = n.burn, cachesize = cachesize, quiet = quiet,
                     param.m = param.m, param.c = param.c, param.a = param.a,
                     param.b = param.b, ...)

  ans <- PeriodCPT.main(ans)
  return(ans)
}

data_value_check.norm <- function(object){
  if(!is.numeric(data.set(object)))
    stop(paste0("Data is invalid for '",distribution(object),"' sampling distribution."))
  if("levels" %in% names(attributes(data.set(object))))
    stop(paste0("Data is invalid for '",distribution(object),"' sampling distribution."))
  if(is.matrix(data.set(object))){
    if(ncol(data.set(object)) > 1){
      stop(paste0("Data is invalid for '",distribution(object),"' sampling distribution."))
    }else{
      data.set(object) <- data.set(object)[,1]
    }
  }
  return(object)
}

