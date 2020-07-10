#set.seed(1)
#x <- ts(sample(c(0,1), size = 240, replace = TRUE), frequency = 24)
#object <- PeriodCPT(data = x, distribution = "bern", n.iter = 100,
#                    quiet = TRUE, inits = "ends", init = "ends")
#probs <- seq(0,1,0.25)


quantile <- function(x, ...){
  quants <- eval(paste0(text = "quantile.",distribution(x),
                   "(object = x, probs = probs)"))
  return(quants)
}


quantile_input_check <- function(object, probs){
  if(class(object) != "pcpt")
    stop("Unexpected class of `object`.")
  if(length(object@date) != 1)
    stop("Argument `object` does not appear to be an output from a PeriodCPT funciton.")
  if(length(summarised(object)) != 1 | any(!summarised(object))){
    object <- summarise_chains(object, all = TRUE)
  }

  if(!is.numeric(probs) | anyNA(probs))
    stop("Probabilities must be a numerical vector between 0 and 1.")
  if(any(probs < 0) | any(probs > 1))
    stop("Probabilities must be a numerical vector between 0 and 1.")

  return(object)
}

segID <- function(object){
  x <- result(object, 1)
  IDs <- matrix(1, nrow = nrow(x),
                ncol = periodlength(object))
  for(i in 1:nrow(IDs)){
    if(x[i, "m"] == 1) next
    for(j in 2:x[i, "m"]){
      st <- x[i,paste0("tau",j-1)] + 1
      ed <- x[i,paste0("tau",j)]
      IDs[i, st:ed] <- j
    }
  }
  return(IDs)
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
attributes(Q.bern.fn) <- list(range = cbind(lower = 0, upper = 1))


Q.pois.fn <- function(x, SSinfo, prob, index = 1, param.prior = NULL){
  if(length(x)>1){
    FN <- rep(NA,length(x))
    for(i in 1:length(x)){
      FN[i] <- Q.pois.fn(x=x[i], SSinfo=SSinfo, prob=prob, index=index,
                         param.prior=param.prior)
    }
    return(FN)
  }
  FN <- NA
  if(index == 1){
    FN <- sum(pgamma(x, shape = SSinfo[,"A"], rate = SSinfo[,"B"]) *
                SSinfo[,"freq"]) - prob * sum(SSinfo[,"freq"])
  }
  return(FN)
}
attributes(Q.pois.fn) <- list(range = cbind(lower = 0, upper = Inf))


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
                    sd = sqrt(SSinfo[,"C"]*param.prior$const.var)) *
                SSinfo[,"freq"]) - prob * sum(SSinfo[,"freq"])
  }
  return(FN)
}
attributes(Q.mean.fn) <- list(range = cbind(lower = -Inf, upper = Inf))


Q.var.fn <- function(x, SSinfo, prob, index = 1, param.prior = NULL){
  if(length(x)>1){
    FN <- rep(NA,length(x))
    for(i in 1:length(x)){
      FN[i] <- Q.var.fn(x=x[i], SSinfo=SSinfo, prob=prob, index=index,
                        param.prior=param.prior)
    }
    return(FN)
  }

  FN <- NA
  if(index == 1){
    FN <- sum(pgamma(1 / x, shape = SSinfo[,"A"], rate = SSinfo[,"B"],lower.tail = FALSE) *
                SSinfo[,"freq"]) - prob * sum(SSinfo[,"freq"])
  }
  return(FN)
}
attributes(Q.var.fn) <- list(range = cbind(lower = 0, upper = Inf))


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
attributes(Q.norm.fn) <- list(range = cbind(lower = c(-Inf,0), upper = c(Inf,Inf)))






quantile.master <- function(object, probs = seq(0, 1, 0.25)){
  object <- quantile_input_check(object, probs)

  quants <- array(NA, dim = c(periodlength(object), length(probs),
                              nsegparam(object)))
  IDs <- segID(object)

  FNs  <- get(paste0("Q.",distribution(object),".fn"))

  SSnames <- names(param.prior(object))
  L <- grepl("param.", SSnames)
  nSuffStats <- sum(L)
  SSnames <- toupper(sub("param.","",SSnames[L]))
  x <- result(object, 1)
  for(i in 1:periodlength(object)){
    for(p in 1:length(probs)){
      PLACE <- p + ((1:nsegparam(object))-1)*length(probs)
      if(probs[p] == 0){
        quants[i,PLACE] <- as.numeric(attributes(FNs)$range[,1])
      }else if(probs[p] == 1){
        quants[i,PLACE] <- as.numeric(attributes(FNs)$range[,2])
      }else{
        SSinfo <- matrix(NA, nrow = nrow(IDs), ncol = 1+nSuffStats)
        for(k in 1:nrow(IDs)){
          SSinfo[k,] <- x[k, c("freq", paste0(SSnames,IDs[k,i]))]
        }
        colnames(SSinfo) <- c("freq", SSnames)

        for(j in 1:nsegparam(object)){
          lim <- attributes(FNs)$range[j,]
          if(lim[1] == -Inf) lim[1] <- -.Machine$double.xmax
          if(lim[2] == +Inf) lim[1] <- +.Machine$double.xmax
          root <- uniroot(f = FNs, SSinfo = SSinfo, prob = probs[p],
                        interval = as.numeric(),
                        param.prior = param.prior(object),
                        index = j)
          quants[i,PLACE[j]] <- root$root
        }
      }
    }
  }
  probtext <- paste0(round(probs*100,5),"%")
  if(nsegparam(object) == 1){
    colnames(quants) <- probtext
  }else{
    colnames(quants) <- paste0(rep(paste0("param",1:3),
                                   each=length(probs)),"_",probtext)
  }
  rownames(quants) <- as.character(1:periodlength(object))
  return(quants)
}

segID <- function(object){
  x <- result(object, 1)
  IDs <- matrix(1, nrow = nrow(x),
                ncol = periodlength(object))
  for(i in 1:nrow(IDs)){
    if(x[i, "m"] == 1) next
    for(j in 2:x[i, "m"]){
      st <- x[i,paste0("tau",j-1)] + 1
      ed <- x[i,paste0("tau",j)]
      IDs[i, st:ed] <- j
    }
  }
  return(IDs)
}







