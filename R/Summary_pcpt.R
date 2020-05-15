Summarise_EvalSS <- function(object, index){
   #index missing -> combine all chains
   #index = 0 -> summarise chains individually
   #index = chainID -> summarise just single chain

  if(missing(index)){
    draws <- NULL
    for(i in 1:n.chains(object)){
      draws <- rbind(draws, MCMC.chain(object, i))
    }
  }else if(index == 0){
     out <- lapply(1:n.chains(object), Summarise_EvalSS, object = object)
     names(out) <- names(MCMC.chains(object))
     return(out)
  }else{
     draws <- MCMC.chain(object, index)
  }

  for(i in ncol(draws):1) draws <- draws[order(draws[,i]),]
  m <- apply(!is.na(draws),1,sum)
  L1 <- c(m[-1]!=m[-length(m)],TRUE)
  L2 <- c(!apply(draws[-1,]==draws[-nrow(draws),],1,all,na.rm=TRUE),TRUE)
  LL <- L1 | L2
  freq <- diff(c(0,which(LL)))
  m <- m[LL]
  draws <- cbind(m,freq,draws[LL,])[order(freq,decreasing = TRUE),]

  phyp_names <- strsplit(names(param.prior(object)),"[.]")
  L <- unlist(lapply(phyp_names,function(a){a[1]=="param"}))
  nSuffStats <- sum(L)
  phyp_names <- unlist(lapply(phyp_names,tail,n=1))

  BLANK <- -1
  offset <- start(data.set(object))[2]
  time <- seq(from = offset, by = 1, length.out = length(data.set(object)))
  time2C <- ((time - 1) %% periodlength(object)) + 1
  drawtmp <- as.numeric(t(draws))
  drawtmp[is.na(drawtmp)] <- BLANK


  C_SuffStats <- .C("PeriodCPT_EvaluateSufficientStats",
        data    = as.numeric(data.set(object)),
        time    = as.integer(time2C),
        n       = as.integer(length(data.set(object))),
        N       = as.integer(periodlength(object)),
        Mdist   = as.character(pcpt.prior(object)$Mprior),
        Pdist   = as.character(distribution(object)),
        Phyp    = as.numeric(param.prior(object)),
        draws   = as.integer(drawtmp),
        ncdraws = as.integer(ncol(draws)),
        nrdraws = as.integer(nrow(draws)),
        nSuffStats = as.integer(nSuffStats),
        blank   = as.numeric(BLANK),
        error   = as.integer(0L),
        out     = vector("numeric", nrow(draws)*nSuffStats*(ncol(draws)-2))
  )

  SufStat <- matrix(C_SuffStats$out,nr=nrow(draws),byrow=TRUE)
  SufStat[SufStat == BLANK] <- NA
  colnames(SufStat) <- paste0(rep(phyp_names,each=ncol(draws)-2),1:(ncol(draws)-2))
  draws <- cbind(draws,SufStat)
  attributes(draws)$info <- list(periodlength = periodlength(object), minseglen = minseglen(object),
                        distribution = distribution(object))
  return(draws)
}

table_npcpt <- function(x){
   N <- attributes(x)$info$periodlength
   l <- attributes(x)$info$minseglen
   maxM <- floor(N/l)
   freq <- rep(0, maxM)
   for(i in 1:nrow(x)){
      freq[x[i,1]] <- freq[x[i,1]] + x[i,2]
   }
   names(freq) <- paste0("m=",1:maxM)
   return(as.table(freq))
}

table_pcpt <- function(x){
   N <- attributes(x)$info$periodlength
   freq <- rep(0,maxM)
   for(i in 1:nrow(x)){
      for(j in 1:x[i,1]){
         freq[x[i,2+j]] <- freq[x[i,2+j]] + x[i,2]
      }
   }
   names(freq) <- paste0("tau",1:N)
   return(as.table(freq))
}

pcpt_mode <- function(x){
  if(is.list(x)){
    out <- lapply(x, pcpt_mode)
    names(out) <- names(x)
    return(x)
  }else{
    mode <- x[1,grepl("tau",colnames(x))]
    return(mode[!is.na(mode)])
  }
}

param_mode <- function(x){
   if(is.list(x)){
      out <- lapply(x, param_mode)
   }else{
      dist <- attributes(x)$info$distribution
      ntau <- sum(grepl("tau",colnames(x)))
      mode_SS <- matrix(x[1,-(1:(2+ntau))],nr=ntau)
      mode_SS <- mode_SS[1:x[1,1], , drop=FALSE]
      for(i in 1:nrow(mode_SS)){
        outi <- eval(parse(text = paste0("param_mode_calc.",
                                      dist,"(mode_SS[i,])")))
        if(i==1){
           out <- matrix(NA,nr=length(outi),nc=nrow(mode_SS))
           rownames(out) <- names(outi)
           colnames(out) <- paste0("seg",1:nrow(mode_SS))
        }
        out[,i] <- outi
      }
   }
   return(out)
}


Evaluate_pcpt_mode <- function(object, summary, combine = TRUE){
   if(combine){
      summary <- Summarise_EvalSS(object)
   }else{
      summary <- Summarise_EvalSS(object, index = 0)
   }

   value <- pcpt_mode(summary)
   if(!is.list(value)){
      value <- as.list(value)
      names(value) <- 1
   }
   pcpt.mode(object) <- value
   return(object)
}

Evaluate_param_mode <- function(object, summary, combine = TRUE){
   if(combine){
      summary <- Summarise_EvalSS(object)
   }else{
      summary <- Summarise_EvalSS(object, index = 0)
   }

   value <- param_mode(summary)
   if(!is.list(value)){
      value <- as.list(value)
      names(value) <- 1
   }
   param.mode(object) <- value
   return(object)
}


