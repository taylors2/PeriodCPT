param.prior.make.norm <- function(param.m = 0, param.c = 1, param.a = 1, param.b = 1, ...){
  ##y ~ Norm(theta[1], theta[2]),
  ##theta[1] ~ Norm(param.m, param.c * theta[2])
  ##theta[2] ~ IGamma(param.a, param.b)

  if(!is.numeric(param.m) || length(param.m) != 1 || anyNA(param.m))
    stop("Hyper-parameter `param.m` specified incorrectly.")
  if(!is.numeric(param.c) || length(param.c) != 1 || any(param.c <= 0) || anyNA(param.c))
    stop("Hyper-parameter `param.c` specified incorrectly.")
  if(!is.numeric(param.a) || length(param.a) != 1 || any(param.a <= 0) || anyNA(param.a))
    stop("Hyper-parameter `param.a` specified incorrectly.")
  if(!is.numeric(param.b) || length(param.b) != 1 || any(param.b <= 0) || anyNA(param.b))
    stop("Hyper-parameter `param.b` specified incorrectly.")
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
  if(missing(data)) stop("Data is missing.")
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

get_Q.norm.fn <- function(){
  return(Q.norm.fn)
}

Q.norm.fn <- function(x, SSinfo, prob, index = 1, param.prior = NULL){

  if(length(x)>1){
    FN <- rep(NA,length(x))
    for(i in 1:length(x)){
      FN[i] <- Q.norm.fn(x=x[i], SSinfo=SSinfo, prob=prob, index=index,
                         param.prior=param.prior)
    }
    return(FN)
  }

  FN <- NA
  if(index == 1){
    y <- (x - SSinfo[,"M"])/sqrt( SSinfo[,"B"]*SSinfo[,"C"]/SSinfo[,"A"]  )
    FN <- sum(pt(y, df = 2*SSinfo[,"A"]) * SSinfo[,"freq"]) -
      prob * sum(SSinfo[,"freq"])
  }else if(index ==2 ){
    FN <- sum(pgamma(1 / x, shape = SSinfo[,"A"], rate = SSinfo[,"B"],lower.tail = FALSE) *
                SSinfo[,"freq"]) - prob * sum(SSinfo[,"freq"])
  }
  return(FN)
}
Q.norm.range <- function(){return(cbind(lower = c(-Inf,0), upper = c(Inf,Inf)))}


