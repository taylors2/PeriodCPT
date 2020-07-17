quantile.master <- function(object, probs = seq(0, 1, 0.25), ...){

  object <- quantile_input_check(object = object, probs = probs)

  quants <- array(NA, dim = c(periodlength(object), length(probs) * nsegparam(object)))
  IDs <- segID(object)

  FNs  <- eval(parse(text = paste0("get_Q.",distribution(object),".fn()")))
  RANGE <- eval(parse(text = paste0("Q.",distribution(object),".range()")))

  SSnames <- names(param.prior(object))
  L <- grepl("param.", SSnames)
  nSuffStats <- sum(L)
  SSnames <- toupper(sub("param.","",SSnames[L]))
  x <- result(object, 1)
  for(i in 1:periodlength(object)){
    for(p in 1:length(probs)){
      PLACE <- p + ((1:nsegparam(object))-1)*length(probs)
      if(probs[p] == 0){
        quants[i,PLACE] <- RANGE[,1]
      }else if(probs[p] == 1){
        quants[i,PLACE] <- RANGE[,2]
      }else{
        SSinfo <- matrix(NA, nrow = nrow(IDs), ncol = 1+nSuffStats)
        for(k in 1:nrow(IDs)){
          SSinfo[k,] <- x[k, c("freq", paste0(SSnames,IDs[k,i]))]
        }
        colnames(SSinfo) <- c("freq", SSnames)

        for(j in 1:nsegparam(object)){
          lim <- RANGE[j,]
          if(lim[1] == -Inf) lim[1] <- -.Machine$double.xmax
          if(lim[2] == +Inf) lim[2] <- +.Machine$double.xmax
          root <- uniroot(f = FNs, SSinfo = SSinfo, prob = probs[p],
                          interval = lim,
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

quantile_input_check <- function(object, probs){
  if(length(object@date) != 1)
    stop("Argument `object` does not appear to be an output from a PeriodCPT function.")
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







