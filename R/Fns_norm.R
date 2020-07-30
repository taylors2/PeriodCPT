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



FNs.norm <- function(x, prob, SSinfo, param.prior = NULL, index = 1){
  q <- rep(NA,length(x))
  for(i in 1:length(x)){
    if(index == 1){
      xx <- (x[i] - SSinfo[,"M"])/sqrt( SSinfo[,"B"]*SSinfo[,"C"]/SSinfo[,"A"]  )
      px <- pt(xx, df = 2*SSinfo[,"A"])
    }else{
      px <- pgamma(1 / x[i], shape = SSinfo[,"A"], rate = SSinfo[,"B"], lower.tail = FALSE)
    }
    p  <- sum(px*SSinfo[,"freq"])/sum(SSinfo[,"freq"])
    q[i] <- p - prob
  }
  return(q)
}
get_FNs.norm <- function(){return(FNs.norm)}
RNG.norm <- function(prob, SSinfo, param.prior = NULL, index = 1){
  if(index == 1){
    q <- qt(prob, df = 2 * SSinfo[,"A"])
    q <- SSinfo[,"M"] + q * sqrt( SSinfo[,"B"]*SSinfo[,"C"]/SSinfo[,"A"] )
  }else{
    q <- 1/qgamma(prob, shape = SSinfo[,"A"], rate = SSinfo[,"B"], lower.tail = FALSE)
  }
  return(range(q))
}

plot.est.format.norm <- function(object, probs = 0.5, param = 0){
  if(param == 1L | param == 2L){
    Q <- quantile(object, probs = probs)
    out <- Q[, grepl(paste0("param",param),colnames(Q)), drop = FALSE]
  }else{
    Q <- quantile(object, probs = 0.5)
    MU <- Q[,1]
    SD <- sqrt(Q[,2])
    out <- matrix(NA, nrow = periodlength(object), ncol = 3)
    out[,1] <- qnorm(probs[1]/2, mean = MU, sd = SD, lower.tail = TRUE )
    out[,2] <- qnorm(0.5,        mean = MU, sd = SD, lower.tail = TRUE )
    out[,3] <- qnorm(probs[1]/2, mean = MU, sd = SD, lower.tail = FALSE)
  }
  return(out)
}

