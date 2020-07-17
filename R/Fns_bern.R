param.prior.make.bern <- function(param.a = 1, param.b = 1, ...){
  ##y~bern(theta), theta~beta(param.a, param.b)
  if(!is.numeric(param.a) || length(param.a) != 1 || any(param.a <= 0) || anyNA(param.a))
    stop("Hyper-parameter `param.a` specified incorrectly.")
  if(!is.numeric(param.b) || length(param.b) != 1 || any(param.b <= 0) || anyNA(param.b))
    stop("Hyper-parameter `param.b` specified incorrectly.")
  out <- c("param.a" = param.a, "param.b" = param.b)
  return(out)
}


PeriodCPT.bern <- function(data, periodlength = NULL, minseglen = 1, Mprior = c("pois", "unif"),
                           Mhyp = 1, spread = 1, param.a = 1, param.b = 1,
                           inits = NULL, n.iter = 1e6, n.chains = 1, n.burn = 0,
                           cachesize = 50, quiet = FALSE, ...){

  distribution <- "bern"
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

data_value_check.bern <- function(object){
  if(!is.numeric(data.set(object)))
    stop(paste0("Data is invalid for '",distribution(object),"' sampling distribution."))
  if("levels" %in% names(attributes(data.set(object))))
    stop(paste0("Data is invalid for '",distribution(object),"' sampling distribution."))
  if(!all(data.set(object) == 0L | data.set(object) == 1L))
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

get_Q.bern.fn <- function(){
  return(Q.bern.fn)
}

Q.bern.fn <- function(x, SSinfo, prob, index = 1, param.prior = NULL){
  if(length(x)>1){
    FN <- rep(NA,length(x))
    for(i in 1:length(x)){
      FN[i] <- Q.bern.fn(x=x[i], SSinfo=SSinfo, prob=prob, index=index,
                         param.prior=param.prior)
    }
    return(FN)
  }
  FN <- NA
  if(index == 1){
    FN <- sum(pbeta(x, shape1 = SSinfo[,"A"], shape2 = SSinfo[,"B"]) *
                SSinfo[,"freq"]) - prob * sum(SSinfo[,"freq"])
  }
  return(FN)
}
Q.bern.range <- function(){return(cbind(lower = 0, upper = 1))}

