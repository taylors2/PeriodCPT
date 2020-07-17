param.prior.make.mean <- function(param.m = 0, param.c = 1, const.var = 1, ...){
  ##y ~ Norm(theta, const.var), theta~Norm(param.m, param.c)
  if(!is.numeric(param.m) || length(param.m) != 1 || anyNA(param.m))
    stop("Hyper-parameter `param.m` specified incorrectly.")
  if(!is.numeric(param.c) || length(param.c) != 1 || any(param.c <= 0) || anyNA(param.c))
    stop("Hyper-parameter `param.c` specified incorrectly.")
  if(!is.numeric(const.var) || length(const.var) != 1 || any(const.var <= 0) || anyNA(const.var))
    stop("Known constant `const.var` is specified incorrectly.")
  out <- c("param.m" = param.m, "param.c" = param.c, "const.var" = const.var)
  return(out)
}

PeriodCPT.mean <- function(data, periodlength = NULL, minseglen = 1, Mprior = c("pois", "unif"),
                           Mhyp = 1, spread = 1, param.m = 0, param.c = 1, const.var = 1,
                           inits = NULL, n.iter = 1e6, n.chains = 1, n.burn = 0,
                           cachesize=50, quiet=FALSE, ...){

  distribution <- "mean"
  nsegparam <- 1
  Mprior <- Mprior[1]
  if(missing(data)) stop("Data is missing.")
  ans <- class_input(data = data, periodlength = periodlength, minseglen = minseglen,
                     distribution = distribution, nsegparam = nsegparam, Mprior = Mprior, Mhyp = Mhyp,
                     spread = spread, inits = inits, n.iter = n.iter,
                     n.chains = n.chains, n.burn = n.burn, cachesize = cachesize,
                     quiet = quiet, param.m = param.m, param.c = param.c,
                     const.var = const.var, ...)

  ans <- PeriodCPT.main(ans)
  return(ans)
}

data_value_check.mean <- function(object){
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

get_Q.mean.fn <- function(){
  return(Q.mean.fn)
}

Q.mean.fn <- function(x, SSinfo, prob, index = 1, param.prior = NULL){
  if(length(x)>1){
    FN <- rep(NA,length(x))
    for(i in 1:length(x)){
      FN[i] <- Q.mean.fn(x=x[i], SSinfo=SSinfo, prob=prob, index=index,
                         param.prior=param.prior)
    }
    return(FN)
  }
  FN <- NA
  if(index == 1){
    FN <- sum(pnorm(x, mean = SSinfo[,"M"],
                    sd = sqrt(SSinfo[,"C"]*param.prior["const.var"])) *
                SSinfo[,"freq"]) - prob * sum(SSinfo[,"freq"])
  }
  return(FN)
}
Q.mean.range <- function(){ return(cbind(lower = -Inf, upper = Inf))}


