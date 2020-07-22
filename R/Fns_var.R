param.prior.make.var <- function(param.a = 1, param.b = 1, ...){
  ##y ~ Norm(0, theta), theta ~ InvGama(param.a, param.b)

  if(!is.numeric(param.a) || length(param.a) != 1 || any(param.a <= 0) || anyNA(param.a))
    stop("Hyper-parameter `param.a` specified incorrectly.")
  if(!is.numeric(param.b) || length(param.b) != 1 || any(param.b <= 0) || anyNA(param.b))
    stop("Hyper-parameter `param.b` specified incorrectly.")
  out <- c("param.a" = param.a, "param.b" = param.b)
  return(out)
}


PeriodCPT.var <- function(data, periodlength = NULL, minseglen = 1, Mprior = c("pois", "unif"),
                          Mhyp = 1, spread = 1, param.a = 1, param.b = 1,
                          inits = NULL, n.iter = 1e6, n.chains = 1, n.burn = 0,
                          cachesize=50, quiet=FALSE, ...){

  distribution <- "var"
  nsegparam <- 1
  Mprior <- Mprior[1]
  if(missing(data)) stop("Data is missing.")
  ans <- class_input(data = data, periodlength = periodlength, minseglen = minseglen,
                     distribution = distribution, nsegparam = nsegparam, Mprior = Mprior, Mhyp = Mhyp,
                     spread = spread, inits = inits, n.iter = n.iter, n.chains = n.chains,
                     n.burn = n.burn, cachesize = cachesize, quiet = quiet,
                     param.a = param.a, param.b = param.b, ...)

  ans <- PeriodCPT.main(ans)
  return(ans)
}

data_value_check.var <- function(object){
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

FNs.var <- function(x, prob, SSinfo, param.prior = NULL, index = 1){
  q <- rep(NA,length(x))
  for(i in 1:length(x)){
    px <- pgamma(1 / x, shape = SSinfo[,"A"], rate = SSinfo[,"B"],lower.tail = FALSE)
    p  <- sum(px*SSinfo[,"freq"])/sum(SSinfo[,"freq"])
    q[i] <- p - prob
  }
  return(q)
}
get_FNs.var <- function(){return(FNs.var)}
RNG.var <- function(prob, SSinfo, param.prior = NULL, index = 1){
  q <- 1 / qgamma(prob, shape = SSinfo[,"A"], rate = SSinfo[,"B"], lower.tail = FALSE)
  return(range(q))
}
