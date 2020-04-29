#############################
##  Wait-time calculation  ##
#############################

WaitTime <- function(DRAW, start=NA, end=NA, prob=0.95, format=c("steps","dys","hrs","mins"), clip=FALSE){
  
  format = match.arg(format)
  N <- attributes(DRAW)$N
  time_txt <- attributes(DRAW)$readtime
  if(!is.na(start) & !is.na(end)){
    if(is.numeric(start)){
      if(!(start %in% (0:N))) stop("Invalid start")  
    }else if(is.character(start)){
      i <- which(time_txt==start)[1]
      if(length(i)==0) stop("Invalid start")
      start <- i
    }else{
      stop("Invalid start")
    }
    
    if(is.numeric(end)){
      if(!(end %in% (0:N))) stop("Invalid end")  
    }else if(is.character(end)){
      i <- which(time_txt==end)[1]
      if(length(i)==0) stop("Invalid end")
      end <- i
    }else{
      stop("Invalid end")
    }
  }else{
    start <- end <- 1
  }
  
  if(start == end){
    start <- 1
    end <- N
  }else{
    end <- (end-2)%%N + 1 ##
  }
  
  waittime = 0;
  if (start<end){
    offsets <- start:end
  }else{
    offsets <- c(start:N,1:end)
  }
  
  for(i in offsets){
    p <- unname(as.numeric(Pred_Prob_of_Next_Obs(DRAW,i)))
    if(clip){
      h <- min(which(p >= prob)[1], length(offsets[-(1:which(offsets==i))]))
    }else{
      h <- which(p>=prob)[1]
    }
    waittime <- max(waittime,h)
  }
  
  wd <- 24*60/N
  if(format == "steps"){
    return(waittime)
  }else if(format == "dys"){
    return(waittime*wd/(24*60))
  }else if(format == "hrs"){
    return(waittime*wd/60)
  }else{# (fromat == "mins")
    return(waittime*wd)
  }
}

WaitTimeTable <- function(DRAW, prob = 0.95){
  N <- attributes(DRAW)$N
  time_txt <- attributes(DRAW)$readtime
  H_step <- H_text <- rep(NA,N)
  MAT <- matrix(NA,N,N)
  for(i in 1:N){
    p <- Pred_Prob_of_Next_Obs(DRAW, offset = i)
    MAT[,i] <- as.numeric(p)
    h <- which(p>=prob)[1]
    if(is.na(h)){
      H_step[i] <- N
      H_text[i] <- time_txt[i]
    }else{
      H_step[i] <- as.numeric(h)
      H_text[i] <- names(h)
    }
  }
  DF.st <- DF.ed <- DF.h <- NULL
  for(t in 1:N){
    i <- which(H_text == time_txt[t])
    if(length(i)>0){
      j <- which.max(H_step[i])
      DF.st <- c(DF.st,time_txt[i[j]])
      DF.ed <- c(DF.ed,H_text[i[j]])
      DF.h  <- c(DF.h, H_step[i[j]])
    }
  }
  out <- data.frame(start = DF.st, end=DF.ed, length=DF.h)
  return(out)
}


