summarize_chains <- function(object, all = TRUE){
  return(summarise_chains(object, all))
}

summarise_chains <- function(object, all = TRUE){

  if(class(object) != "pcpt") stop("Unexpected class of `object`.")
  if(length(object@date) == 0)
    stop("Argument `object` does not appear to be an output from a PeriodCPT function.")
  if(length(summarised(object)) == 1 & all(summarised(object)))
    return(object)  #Already summarised & combinded!!!
  if(length(all)!=1 | !is.logical(all) | anyNA(all))
    stop("Argument `all` is not a single logical value.")

  for(index in 1:n.chains(object)){ ##Perform summary per chain
    object <- summarise_single_chain(object, index)
  }
  if(all) object <- summarise_combine(object)
  return(object)
}

#summarize_combine <- function(object){
#  return(summarise_combine(object))
#}

summarise_combine <- function(object){

  #if(any(!summarised(object))){
  #  ##There are some chains that are not summarised?
  #  return(summarise_chains(object, all = TRUE))
  #}
  if(n.chains(object) == 1) #There is only one summarised chain, done!
    return(object)
#  if(length(summarised(object)) == 1 & all(summarised(object))) #Summarised chains have already been combined
#    return(object)

  tab <- result(object, 1)
  for(c in 2:n.chains(object)){
    tab2 <- result(object, 2)
    ##Pad tables such that ncol(tab)==ncol(tab2)
    if(ncol(tab) < ncol(tab2)){
      tab_tmp <- array(NA,dim=c(nrow(tab), ncol(tab2)))
      tab_tmp[,colnames(tab2) %in% colnames(tab)] <- tab
      colnames(tab_tmp) <- colnames(tab2)
      tab <- tab_tmp
    }else if(ncol(tab) > ncol(tab2)){
      tab_tmp <- array(NA,dim=c(nrow(tab2), ncol(tab)))
      tab_tmp[,colnames(tab) %in% colnames(tab2)] <- tab2
      colnames(tab_tmp) <- colnames(tab)
      tab2 <- tab_tmp
    }

    keep <- rep(TRUE,nrow(tab2)) #tab2 rows not in tab == TRUE (needed for appending)
    for(i in 1:nrow(tab2)){
      j <- 1
      while(j<=nrow(tab) && keep[i]){
        if(tab2[i,1]!=tab[j,1]){ #number of pcpts dont match
          j <- j+1
        }else if(all(tab2[i,2+tab2[i,1]] == tab[j,2+tab[j,1]])){
          #pcpts match, add freq to tab and mark as combined
          tab[j,2] <- tab[j,2] + tab2[i,2] ##Add frequencies
          keep[i] <- FALSE
        }else{
          j <- j+1
        }
      }
    }
    ##Append uniqe draw from tab2 into tab.
    tab <- rbind(tab, tab2[keep,,drop=FALSE])
  }

  tab <- tab[order(tab[,2],decreasing = TRUE),, drop=FALSE]

  #Put tab matrix into required format for results slot and assign
  out <- list(tab)
  names(out) <- 1
  results(object) <- out
  summarised(object) <- TRUE
  return(object)
}

#summarize_single_chain <- function(object, index){
#  return(summarise_single_chain(object, index))
#}

summarise_single_chain <- function(object, index = 1){

#  if(length(index)!=1 | !is.numeric(index) | anyNA(index))
#    stop("Argument `index` is not a single integer value for a identifying chain.")
#  if(index<=0 | index>n.chains(object) | floor(index)!=index)
#    stop("Argument `index` is not a single integer value for a identifying chain.")

  if(summarised(object)[index]) return(object) #chain already summarised!!!


  ##Tabulating chain samples
  chain <- result(object, index)
  ntaumax <- ncol(chain)
  niter <- nrow(chain)
  chain <- as.numeric(t(chain))
  BLANK <- -1
  chain[is.na(chain)] <- BLANK
  tally <- .C("Tally_pcpt", as.integer(chain), as.integer(niter), as.integer(ntaumax),
              as.integer(BLANK), tally = vector("integer", length = niter))$tally
  uniquedraws <- result(object,index)[tally>0,,drop=FALSE]
  tally <- tally[tally>0]
  m <- apply(!is.na(uniquedraws),1,sum)

  ##Evaluating Sufficient Stats
  SSnames <- names(param.prior(object))
  L <- grepl("param.", SSnames)
  nSuffStats <- sum(L)
  phyp_names <- toupper(sub("param.","",SSnames[L]))

  offset <- start(data.set(object))[2]
  time <- seq(from = offset, by = 1, length.out = length(data.set(object)))
  time2C <- ((time - 1) %% periodlength(object)) + 1
  drawtmp <- as.numeric(t(uniquedraws))
  drawtmp[is.na(drawtmp)] <- BLANK

  C_SuffStats <- .C("PeriodCPT_EvaluateSufficientStats",
                    data    = as.numeric(data.set(object)),
                    time    = as.integer(time2C),
                    n       = as.integer(length(data.set(object))),
                    N       = as.integer(periodlength(object)),
                    Pdist   = as.character(distribution(object)),
                    Phyp    = as.numeric(param.prior(object)),
                    draws   = as.integer(drawtmp),
                    ncdraws = as.integer(ncol(uniquedraws)),
                    nrdraws = as.integer(nrow(uniquedraws)),
                    nSuffStats = as.integer(nSuffStats),
                    blank   = as.integer(BLANK),
                    error   = as.integer(0L),
                    out     = vector("numeric", nSuffStats*length(uniquedraws))
  )

  SufStat <- matrix(C_SuffStats$out,nrow=nrow(uniquedraws),byrow=TRUE)
  SufStat[SufStat == BLANK] <- NA
  tab <- cbind(m,tally,uniquedraws,SufStat)
  colnames(tab) <- c("m", "freq", paste0(rep(c("tau",phyp_names), each=ncol(uniquedraws)),
                                         1:ncol(uniquedraws)))

  tab <- tab[order(tab[,2],decreasing = TRUE),, drop=FALSE]

  ####Assign summary table back to appropriate result slot
  result(object, index) <- tab
  summarise_logical <- summarised(object)
  summarise_logical[index] <- TRUE
  summarised(object) <- summarise_logical

  return(object)
}

table_npcpt <- function(object){
  if(class(object) != "pcpt") stop("Unexpected class of `object`.")
  if(length(object@date) == 0)
    stop("Argument `object` does not appear to be an output from a PeriodCPT function.")

  len_result <- length(summarised(object))
  tab <- matrix(0,nrow=len_result, ncol=npcpts.max(object))
  for(index in 1:len_result){
    if(summarised(object)[index]){
      x <- result(object, index)
      for(m in 1:npcpts.max(object)){
        tab[index,m] <- sum(x[x[, "m"]==m, "freq"])
      }
    }else{
      m_chain <- apply(!is.na(result(object, index)),1,sum)
      tab[index,] <- as.numeric(table(c(m_chain,1:ncol(tab))))-1
    }
  }

  if(len_result > 1){
    colnames(tab) <- paste0("m",1:ncol(tab))
    rownames(tab) <- paste0("Chain ",1:nrow(tab))
  }else{
    tab <- as.numeric(tab)
    names(tab) <- paste0("m",1:length(tab))
  }
  return(as.table(tab))
}


table_pcpt <- function(object){

  if(class(object) != "pcpt") stop("Unexpected class of `object`.")
  if(length(object@date) == 0)
    stop("Argument `object` does not appear to be an output from a PeriodCPT function.")
  N <- periodlength(object)

  len_result <- length(summarised(object))
  tab <- matrix(0,nrow=len_result, ncol=N+1)
  for(index in 1:len_result){
    x <- result(object, index)
    if(summarised(object)[index]){
      for(i in 1:nrow(x)){
        if(x[i,"m"] == 1){
          tab[index, N+1] <- x[i, "freq"]
        }else{
          for(j in 1:x[i,"m"]){
            tau <- x[i, paste0("tau",j)]
            tab[index, tau] <- tab[index, tau] + x[i, "freq"]
          }
        }
      }
    }else{
      L <- apply(!is.na(x),1,sum) == 1
      tab[index, N+1] <- sum(L)
      x <- x[!L, , drop=FALSE]
      if(nrow(x)>0){
        for(tau in 1:N){
          tab[index, tau] <- sum(apply(x == tau,1,any,na.rm=TRUE))
        }
      }
    }
  }
  if(len_result > 1){
    colnames(tab) <- c(paste0("tau",1:N),"m1")
    rownames(tab) <- paste0("Chain ",1:nrow(tab))
  }else{
    tab <- as.numeric(tab)
    names(tab) <- c(paste0("tau",1:N),"m1")
  }
  return(as.table(tab))
}



##---------------------------------------------------------------------
quantile_per_time_slot <- function(object, ...){
  cat("Function 'quantile_per_time_slot' not yet implemented\n")
  return(NULL)
}
##---------------------------------------------------------------------
#plot_chain <- function(object, col, plot = TRUE){
#
#
#
#  if(!any(grepl("freq", names(result(object,1))))){
#    stop("Cannot create plot on summarised chain.")
#  }
#  tmp <- floor(log10(n.iter(object)))
#  if(tmp > 2){
#    by = 10^(tmp-2)
#    rng <- seq(from = 10^(tmp-2), to = n.iter(object), by = by)
#  }else{
#    rng <- 1:n.iter(object)
#  }
#  N <- periodlength(object)
#  out <- vector("list", n.chains(object))
#  for(c in 1:n.chains(object)){
#    mat <- matrix(0,nr=length(rng),nc=N+1)
#    colnames(mat) <- c(paste0("tau",1:N),"m1")
#    rownames(mat) <- rng
#    j <- 1
#    mcmc <- result(object, c)
#    for(i in 1:n.iter(object)){
#      tau <- mcmc[i, !is.na(mcmc[i,])]
#      if(length(tau)==1){
#        mat[j,N+1] <- mat[j,N+1] + 1
#      }else{
#        mat[j,tau] <- mat[j,tau] + 1
#      }
#      if(i == rng[j]){
#        j <- j + 1
#      }
#    }
#
#    out[[c]] <- mat/diff(c(0,rng))
#  }
#  names(out) <- 1:n.chains(out)
#
#  if(plot){
#    N <- periodlength(object)
#    n <- n.chains(object)
#
#    tmp3 <- matrix(0,nr=nrow(mat1),nc=n*ncol(mat1))
#    for(i in 1:n){
#      tmp <- out[[i]]
#      tmp3[,((1:ncol(tmp)) - 1) * n + i ] <- tmp  + 2*(i-1)*(tmp!=0)
#    }
#
#    if(missing(col)) col = 1:n
#    if(!is.character(COLS)){
#      COLS <- c("black","red","green","blue","cyan","magenta",
#                "yellow","grey")[(col - 1)%%8 + 1]
#    }else{
#      COLS <- col
#    }
#    colhexchar <- as.character(as.hexmode(col2rgb(COLS)))
#    for(i in 1:length(colhexchar)){
#      if(nchar(colhexchar[i]) == 1)
#        colhexchar[i] <- paste0("0",colhexchar[i])
#    }
#    hexcol <- paste0("#",apply(colhexchar,2,paste0,collapse=""))
#    UseCOLS <- NULL
#    levels <- as.character(as.hexmode(0:15))
#    for(i in 1:length(COLS)){
#      UseCOLS <- c(UseCOLS,paste0(hexcol[i],levels,levels),
#                   paste0(hexcol[i],"ff"),
#                   paste0("#000000",levels[-1],levels[-1]))
#    }
#    collegend <- paste0("#000000",levels,levels)
#
#    AT <- (0:N)*n
#    rng <- seq(from=0,by=100,length.out = nrow(tmp3))
#    fields::image.plot(rng,1:((N+1)*n),tmp3,zlim=c(0,n*2),col=UseCOLS,
#                       axis.args = list(col="white", at = c(0,n*2), labels=c("","")),
#                       xlab = "iteration", yaxt="n",ylab="")
#    axis(at = seq(from = mean(AT[1:2])+0.5, by = AT[2]-AT[1], len = N+1),
#         side=2, labels=c(paste0("tau",1:N),"m1"), las = 1)
#    fields::image.plot(NA,zlim=c(0,1),col=c("#FFFFFF","#FFFFFF"),
#                       legend.only = TRUE, axis.args = list(col="white", at = c(0,1),
#                          labels=c("","")))  ##To blank out legend
#    fields::image.plot(NA,zlim=c(0,1),col=collegend,legend.only = TRUE)
#    segments(diff(range(rng))*-0.005,AT + 0.5,rep(max(rng)*0.827,13),
#             AT + 0.5)
#    invisible(out)
#  }else{
#    return(out)
#  }
#}


