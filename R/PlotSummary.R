#######################################
##  PLOTTING AND SUMMARTY FUNCTIONS  ##
#######################################

#Circ_plot <- function(data, time=NULL, rlim=NULL, tmax=NULL, dot=TRUE, clocktime=TRUE,
#                      main=NULL, tlab="Time", rlab=NULL,cex.lab=1,cex.axis=1,cex.main=1,...){
#
#  ##Check and format input
#  if(missing(data)) stop("data not defined.")
#  if(is.null(rlim)){
#    rlim <- range(data,na.rm=TRUE)
#  }
#  if(!is.numeric(rlim) | length(rlim)!=2) stop("Invalid rlim.")
#  if(is.null(time)){
#    time <- 1:length(data)
#  }else{
#    if(!is.numeric(time)) stop("Invalid integer vector time.")
#    if(length(time)!=length(data) || any(time<=0) || any(time%%1!=0))
#      stop("Invalid integer vector time.")
#  }
#  if(is.null(tmax)) tmax <- max(time)
#  tmax <- max(tmax,time)
#  if(is.null(rlab)) rlab <- deparse(substitute(data))
#
#  if(clocktime){
#    taxis <- 0
#  }else if(tmax%%8==0){
#    taxis <- 8
#  }else if(tmax%%6==0){
#    taxis <- 6
#  }else if(tmax%%4==0){
#    taxis <- 4
#  }else{
#    taxis <- 8
#  }
#  rtick <- pretty(rlim)
#  rlim <- range(rtick)
#  rtick <- (rtick-rlim[1])/(rlim[2]-rlim[1])  #standardize to 0,1 range
#  dr <- rtick[2]-rtick[1]
#  Rmax <- max(rtick) + dot*dr
#
#  ##Create plotting frame
#  plot(NA,NA,xlim=c(-1,1)*(Rmax+1.5*dr),ylim=c(-1,1)*(Rmax+1.5*dr),
#       ylab="",xlab="", xaxt="n", yaxt="n",frame=FALSE,asp=TRUE)
#  title(main,cex=cex.main)
#
#  ##Grid
#  segments(-Rmax,0,Rmax,0,lty=2,col="grey60")
#  if(taxis == 4){
#    segments(0,-Rmax,0,Rmax,lty=2,col="grey60")
#  }else if(taxis == 6){
#    segments(-Rmax*sin(pi/3),-Rmax*cos(pi/3),Rmax*sin(pi/3),Rmax*cos(pi/3),lty=2,col="grey60")
#    segments(-Rmax*sin(pi/3),Rmax*cos(pi/3),Rmax*sin(pi/3),-Rmax*cos(pi/3),lty=2,col="grey60")
#  }else{
#    segments(0,-Rmax,0,Rmax,lty=2,col="grey60")
#    segments(-Rmax/sqrt(2),-Rmax/sqrt(2),Rmax/sqrt(2),Rmax/sqrt(2),lty=2,col="grey60")
#    segments(-Rmax/sqrt(2),Rmax/sqrt(2),Rmax/sqrt(2),-Rmax/sqrt(2),lty=2,col="grey60")
#  }
#
#  psi <- seq(0,2*pi,len=101)
#  for(r in rtick + dot*dr){
#    polygon(r*cos(psi),r*sin(psi),col=NA,lty=2,border="grey60")
#  }
#
#  #Central Dot
#  if(dot) polygon(0.8*dr*cos(psi),0.8*dr*sin(psi),col="grey90")
#
#  ##Axis Values
#  text(rtick + dot*dr,rep(0,length(rtick)),
#    seq(rlim[1],rlim[2],len=length(rtick)),cex=cex.axis)
#  if(clocktime){
#    text(c(0,0,-1)*(Rmax+dr),c(1,-1,0)*(Rmax+dr),
#         paste0(c("00","12","18"),":00"))
#    text(c(1,1,-1,-1)*(Rmax+dr)/sqrt(2),
#         c(1,-1,-1,1)*(Rmax+dr)/sqrt(2),
#         paste0(c("03","09","15","21"),":00"))
#  }else{
#    ang <- pi/2-2*pi*seq(0,1,len=taxis+1)[-1]
#    txt <- paste0(1:taxis,"/",taxis)
#    if(taxis==4){ang <- ang[-1];txt <- txt[-1]}
#    if(taxis==8){ang <- ang[-2];txt <- txt[-2]}
#    text((Rmax+dr)*cos(ang),(Rmax+dr)*sin(ang),txt,cex=cex.axis)
#  }
#
#  ##Axis Names
#  text(0,Rmax+1.5*dr,tlab,cex=cex.lab)
#  text(Rmax,-0.125*Rmax,rlab,pos=2,offset=-1,cex=cex.lab)
#
#  #Plot data
#  R <- (data-rlim[1])/(rlim[2]-rlim[1])
#  R[R<0 | R>Rmax] <- NA  #NULL data outside plotting range
#  R <- R + dot*dr
#  angle <- pi/2 - 2*pi*((time-1)/tmax)
#  points(R*cos(angle), R*sin(angle),...)
#
#  #invisibly return plotting info
#  return(invisible(list(tmax=tmax, Rmax=Rmax, rlim=rlim, dot=dot, rtick=rtick)))
#}
#
#Circ_rline <- function(t, r1=NULL, r2=NULL, plotinfo, ...){
#  if(missing(t)) stop("Missing t.")
#  if(is.character(t)){
#    mins <- sum(c(60,1)*as.numeric(strsplit(t,":")[[1]]))
#    t <- mins*plotinfo$tmax/1440
#    ang <- pi/2 - 2*pi*t/plotinfo$tmax
#  }else{
#    ang <- pi/2 - 2*pi*(t-1)/plotinfo$tmax
#  }
#  if(is.null(r1)){
#    r1 <- rep(0,length(t))
#  }else{
#    r1 <- (r1-plotinfo$rlim[1])/(plotinfo$rlim[2]-plotinfo$rlim[1])
#  }
#  if(is.null(r2)){
#    r2 <- rep(plotinfo$Rmax,length(t))
#  }else{
#    r2 <- (r2-plotinfo$rlim[1])/(plotinfo$rlim[2]-plotinfo$rlim[1])
#  }
#  r1 <- pmin(pmax(r1,0),plotinfo$Rmax)
#  r2 <- pmin(pmax(r2,0),plotinfo$Rmax)
#
#  if(plotinfo$dot){
#    dr <- plotinfo$rtick[2]-plotinfo$rtick[1]
#    r1 <- r1 + dr
#    r2 <- r2 + dr
#  }
#  segments(r1*cos(ang),r1*sin(ang),r2*cos(ang),r2*sin(ang),...)
#}
#
#Circ_arc <- function(r, t1 = NULL, t2 = NULL, plotinfo, ...){
#
#  if(length(t1) > 1 || length(t2) > 1) stop("Invalid t1 and/or t2.")
#  if(missing(r)) stop("Missing r.")
#  if(!is.numeric(r) | length(r)!=1) stop("Invalid r.")
#  r <- (r-plotinfo$rlim[1])/(plotinfo$rlim[2]-plotinfo$rlim[1])
#  r <- pmin(pmax(r,0),plotinfo$Rmax)
#  if(plotinfo$dot) r <- r + (plotinfo$rtick[2]-plotinfo$rtick[1])
#  if(is.null(t1) | is.null(t2)){
#    t1 = 1
#    t2 = plotinfo$tmax+1
#  }
#  if(is.character(t1)){
#    mins1 <- sum(c(60,1)*as.numeric(strsplit(t1,":")[[1]]))
#    t1 <- mins1*plotinfo$tmax/1440 + 1
#  }else if(t1%%1 != 0 || t1<0 || t1>plotinfo$tmax){
#    stop("Invalid t1 and/or t2.")
#  }
#  if(is.character(t2)){
#    mins2 <- sum(c(60,1)*as.numeric(strsplit(t2,":")[[1]]))
#    t2 <- mins2*plotinfo$tmax/1440 + 1
#  }else if(t2%%1 != 0 || t2<0 || t2>(1+plotinfo$tmax)){
#    stop("Invalid t1 and/or t2.")
#  }
#  if(t1<=t2){
#    theta = seq(t1,t2,len=101)-1
#  }else{
#    theta <- NULL
#    if(t1-1 < plotinfo$tmax){
#      theta = c(theta,seq(t1-1,plotinfo$tmax,len=101))
#    }
#    if(t2 > 0) theta <- c(theta,seq(0,t2-1,len=101))
#  }
#  ang = pi/2 - (theta/plotinfo$tmax)*2*pi
#  lines( r*cos(ang), r*sin(ang), ...)
#
#}
#
#Circ_line <- function(radius, time, close=TRUE, plotinfo, ...){
#  radius <- (radius-plotinfo$rlim[1])/(plotinfo$rlim[2]-plotinfo$rlim[1])
#  radius[radius<0 | radius>plotinfo$Rmax] <- NA
#  if(plotinfo$dot) radius <- radius + (plotinfo$rtick[2]-plotinfo$rtick[1])
#  if(missing(time)) time = 1:length(radius)
#  if(is.character(time)){
#    t <- rep(NA,length(time))
#    for(tt in 1:length(time)){
#      mins  <- sum(c(60,1)*as.numeric(strsplit(time[tt],":")[[1]]))
#      t[tt] <- mins*plotinfo$tmax/1440
#    }
#    time <- t+1
#  }
#  if(length(time)%%plotinfo$tmax != 0 && close)
#    warning("line may not be a closed loop.")
#  ang = pi/2 - 2*pi*(time-1)/plotinfo$tmax
#  x = radius*cos(ang)
#  y = radius*sin(ang)
#  lines(x,y,...)
#  if(close){
#    i = c(1,length(x))
#    lines(x[i], y[i],...)
#  }
#}
#
#diagnostic.plot <- function(DRAW, night=NULL, NAME){
#  if(!any("Freq"==colnames(DRAW$TAU))){
#    TAU = SummariseOutput(DRAW)
#  }else{
#    TAU = DRAW$TAU
#  }
#  M = sort(unique(TAU[,1]))
#  Mcount = rep(0,length(M))
#  Mmode = rep(0,length(M))
#  for(i in 1:length(Mcount)){
#    Mcount[i] = sum(TAU[TAU[,1]==M[i],2])
#    Mmode[i] = which(TAU[,1]==M[i])[1]
#  }
#  N <- attributes(DRAW)$N
#  if(missing(NAME)) NAME = attributes(DRAW)$name
#
#  layout(cbind(c(1,1,2),3))
#  plot(NA,NA,ylim=c(0,N),xlim=range(M),xlab="No. of Segments",
#       ylab="Mode Cpt",main=NAME,yaxt="n")
#  axis(side=2,at=(0:4)*(N/4),
#       labels=c("00:00","06:00","12:00","18:00","00:00"))
#  abline(h=(0:4)*(N/4), lty=3, col="grey60")
#  abline(v=M, lty=3, col="grey60")
#
#  for(i in 1:length(Mmode)){
#    if(M[i]>1) points(rep(M[i],M[i]), TAU[Mmode[i],2+(1:M[i])]-1,pch=4)
#  }
#  barplot(Mcount,names.arg=M,ylab="Freq.",xlab="No. of Segments",main="No. of Segments Posterior")
#
#  if(attributes(DRAW)$DIST=="Normal"){
#    E_count <- Param_quantiles(probs=0.5, DRAW=DRAW)
#    tmp <- sqrt(E_count[,2])*qnorm(0.975)
#    E_count <- cbind(E_count[,1],E_count[,1]-tmp,E_count[,1]+tmp)
#  }else{
#    E_count <- Param_quantiles(probs=c(0.5,0.025,0.975), DRAW=DRAW)
#  }
#  if(attributes(DRAW)$DIST=="Var") E_count <- sqrt(E_count)
#  info = Circ_plot(data=DRAW$data, time=DRAW$time, main=NAME)
#  Circ_line(radius=E_count[,1], plotinfo=info, col=2, lwd=2)
#  Circ_line(radius=E_count[,2], plotinfo=info, col=4, lwd=2, lty=2)
#  Circ_line(radius=E_count[,3], plotinfo=info, col=4, lwd=2, lty=2)
#  if(!is.null(night)){
#    if(is.character(night)){
#      wd <- 24*60/N
#      time_txt <- format(seq.POSIXt(from=as.POSIXct("1970-01-01 00:00"),
#                                    to=as.POSIXct("1970-01-02 00:00:00"),
#                                    by = wd*60),"%H:%M")#[-1]
#      time_txt <- time_txt[-length(time_txt)]
#      i1 <- which(time_txt==night[1])[1]
#      i2 <- which(time_txt==night[2])[1]
#      night <- c(i1,i2)
#    }
#    Circ_rline(t=as.numeric(night), plotinfo = info, col=c("darkgreen","orange"),
#               lwd=2, lty=2)
#  }
#}
#
#Param_quantiles = function(probs=0.5, DRAW){
#  if(!is.numeric(probs)) stop("Invalid probabilities.")
#  if(any(probs<=0 | probs>=1)) stop("Invalid probabilities.")
#  if(!any("Freq"==colnames(DRAW$TAU))){
#    TAU = SummariseOutput(DRAW)
#  }else{
#    TAU = DRAW$TAU
#  }
#  if(attributes(DRAW)$DIST=="Poisson"){
#    UPPER <- max(DRAW$data)
#    LOWER <- 0
#    param_dist <- function(x,ABW,p) return(sum(pgamma(x,ABW[,1],ABW[,2])*ABW[,3])-p)
#  }else if(attributes(DRAW)$DIST=="Binomial"){
#    UPPER <- 1
#    LOWER <- 0
#    param_dist <- function(x,ABW,p) return(sum(pbeta(x,ABW[,1],ABW[,2])*ABW[,3])-p)
#  }else if(attributes(DRAW)$DIST=="Normal"){
#    return(Param_quantiles_normal(probs=probs,DRAW=DRAW))
#  }else if(attributes(DRAW)$DIST=="Mean"){
#    UPPER <- max(abs(DRAW$data))
#    LOWER <- -UPPER
#    param_dist <- function(x,ABW,p) return(sum(pnorm(x,mean=ABW[,1],sd=sqrt(ABW[,2]))*ABW[,3])-p)
#  }else if(attributes(DRAW)$DIST=="Var"){
#    UPPER <- max(abs(DRAW$data))
#    LOWER <- 0
#    param_dist <- function(x,ABW,p) return(sum((1-pgamma(1/x,shape=ABW[,1],rate=ABW[,2]))*ABW[,3])-p)
#  }else{
#    stop("Distributional assumption not recognised.")
#  }
#  N = attributes(DRAW)$N
#  E_Quantiles <- matrix(NA,nrow=N,ncol=length(probs))
#  W <- TAU[,2]/sum(TAU[,2])
#  MMax <- max(TAU[,1])
#  for(i in 1:N){
#    A <- B <- rep(NA,nrow(TAU))
#    for(k in 1:nrow(TAU)){
#      if(TAU[k,1]==1){
#        A[k] = TAU[k,2+MMax+1]
#        B[k] = TAU[k,2+2*MMax+1]
#      }else if((i<=TAU[k,3])|(i>TAU[k,2+TAU[k,1]])){
#        A[k] = TAU[k,2+MMax+1]
#        B[k] = TAU[k,2+2*MMax+1]
#      }else{
#        j=2
#        while(!(i<=TAU[k,2+j] & i>TAU[k,1+j])) j = j+1
#        A[k] = TAU[k,2+MMax+j]
#        B[k] = TAU[k,2+2*MMax+j]
#      }
#    }
#    for(p in 1:length(probs)){
#      fUP <- param_dist(x=UPPER,ABW=cbind(A,B,W),p=probs[p])
#      fLO <- param_dist(x=LOWER,ABW=cbind(A,B,W),p=probs[p])
#      while(sign(fLO)!=-sign(fUP)){
#        UPPER <- 2*UPPER
#        LOWER <- 2*LOWER
#        fUP <- param_dist(x=UPPER,ABW=cbind(A,B,W),p=probs[p])
#        fLO <- param_dist(x=LOWER,ABW=cbind(A,B,W),p=probs[p])
#      }
#      u <- uniroot(f=param_dist,interval=c(LOWER,UPPER),
#                   ABW=cbind(A,B,W),p=probs[p])
#      E_Quantiles[i,p] <- u$root
#    }
#  }
#  if(length(probs)==1) return(as.numeric(E_Quantiles))
#  colnames(E_Quantiles) = paste0(100*probs,"%")
#  return(E_Quantiles)
#}
#
#Param_quantiles_normal <- function(probs=0.5,DRAW){
#
#  param_dist_mean <- function(x,MSdW,p) return(sum(pt((x-MSdW[,1])/MSdW[,2],df=MSdW[,3])*MSdW[,4])-p)
#  param_dist_var  <- function(x,ABW,p) return(sum((1-pgamma(1/x,shape=ABW[,1],rate=ABW[,2]))*ABW[,3])-p)
#  if(!any("Freq"==colnames(DRAW$TAU))){
#    TAU = SummariseOutput(DRAW)
#  }else{
#    TAU = DRAW$TAU
#  }
#  N = attributes(DRAW)$N
#  Mean_Quantiles <- Var_Quantiles <- matrix(NA,nrow=N,ncol=length(probs))
#  W <- TAU[,2]/sum(TAU[,2])
#  MMax <- max(TAU[,1])
#  for(i in 1:N){
#    A <- B <- M <- C <- rep(NA,nrow(TAU))
#    for(k in 1:nrow(TAU)){
#      if(TAU[k,1]==1){
#        M[k] = TAU[k,2+MMax+1]
#        C[k] = TAU[k,2+2*MMax+1]
#        A[k] = TAU[k,2+3*MMax+1]
#        B[k] = TAU[k,2+4*MMax+1]
#      }else if((i<=TAU[k,3])|(i>TAU[k,2+TAU[k,1]])){
#        M[k] = TAU[k,2+MMax+1]
#        C[k] = TAU[k,2+2*MMax+1]
#        A[k] = TAU[k,2+3*MMax+1]
#        B[k] = TAU[k,2+4*MMax+1]
#      }else{
#        j=2
#        while(!(i<=TAU[k,2+j] & i>TAU[k,1+j])) j = j+1
#        M[k] = TAU[k,2+MMax+j]
#        C[k] = TAU[k,2+2*MMax+j]
#        A[k] = TAU[k,2+3*MMax+j]
#        B[k] = TAU[k,2+4*MMax+j]
#      }
#    }
#    ABW <- cbind(A,B,W)
#    MSdW <- cbind(M,sqrt(C*B/A),2*A,W)
#    UPPER <- max(abs(DRAW$data))
#    for(p in 1:length(probs)){
#      fLO <- param_dist_mean(-UPPER,MSdW=MSdW,p=probs[p])
#      fUP <- param_dist_mean( UPPER,MSdW=MSdW,p=probs[p])
#      while(sign(fLO)!=-sign(fUP)){
#        UPPER <- UPPER*2
#        fLO <- param_dist_mean(-UPPER,MSdW=MSdW,p=probs[p])
#        fUP <- param_dist_mean( UPPER,MSdW=MSdW,p=probs[p])
#      }
#      u <- uniroot(f=param_dist_mean,interval=c(-1,1)*UPPER,
#                   MSdW=MSdW,p=probs[p])
#      Mean_Quantiles[i,p] <- u$root
#
#      fLO <- param_dist_var(     0,ABW=ABW,p=probs[p])
#      fUP <- param_dist_var( UPPER,ABW=ABW,p=probs[p])
#      while(sign(fLO)!=-sign(fUP)){
#        UPPER <- UPPER*2
#        fUP <- param_dist_var( UPPER,ABW=ABW,p=probs[p])
#      }
#      u <- uniroot(f=param_dist_var,interval=c(0,UPPER),
#                   ABW=ABW,p=probs[p])
#      Var_Quantiles[i,p] <- u$root
#    }
#  }
#  colnames(Mean_Quantiles) <- paste0("mean",100*probs,"%")
#  colnames(Var_Quantiles) <- paste0("var",100*probs,"%")
#  return(cbind(Mean_Quantiles,Var_Quantiles))
#}
#
#CPT_mode <- function(DRAW){
#  if(!any("Freq"==colnames(DRAW$TAU))){
#    TAU = SummariseOutput(DRAW)
#  }else{
#    TAU = DRAW$TAU
#  }
#  if(TAU[1,1]==1) return("NULL")
#  mode = attributes(DRAW)$readtime[TAU[1,2+(1:TAU[1,1])]]
#  return(mode)
#}
#
#CPT_freq <- function(DRAW, marg=TRUE){
#  if(!any("Freq"==colnames(DRAW$TAU))){
#    TAU = SummariseOutput(DRAW)
#  }else{
#    TAU = DRAW$TAU
#  }
#  MMax <- max(TAU[,1])
#  if(marg){
#    N <- attributes(DRAW)$N
#    FREQ <- rep(0,attributes(DRAW)$N+1)
#    for(i in 1:N){
#      tmp <- matrix(TAU[,2+(1:MMax)],ncol=MMax,nrow=nrow(TAU))
#      FREQ[i] = sum(TAU[apply(tmp==i,1,any,na.rm=TRUE),2])
#    }
#    m1 <- which(TAU[,1]==1)
#    if(length(m1)>0){
#      FREQ[1] = FREQ[1] - TAU[m1,2]
#      FREQ[N+1] = TAU[m1,2]
#    }
#    TAB <- as.table(FREQ/sum(TAU[,2]))
#    names(TAB) <- c(attributes(DRAW)$readtime,"NULL")
#  }else{
#    textTAU <- matrix(TAU[,2+(1:MMax)],nrow=nrow(TAU),ncol=MMax)
#    for(i in 1:MMax) textTAU[,i] <- attributes(DRAW)$readtime[TAU[,2+i]]
#    m1 <- which(TAU[,1]==1)
#    if(length(m1)>1) textTAU[m1,1] <- "NULL"
#    textcpt <- apply(textTAU,1,paste0,collapse=", ")
#    TAB <- as.table(TAU[,2]/sum(TAU[,2]))
#    names(TAB) <- textcpt
#  }
#  return(TAB)
#}
#
#
#
#
#QuantileGivenCPT <- function(probs=0.5, cpts = NULL, data = NULL, time = NULL,
#        l = 4, hyp = NULL, DIST = c("Poisson","Binomial","Normal","Mean","Var")){
#
#  DIST <- match.arg(DIST)
#  if(anyNA(data) | anyNA(time)) stop("data and time connot contain missing values.")
#  if(!is.numeric(data) | !is.numeric(time)) stop("data and time must be a numerical vector.")
#  if(length(data)!=length(time)) stop("data and time must have the same length,")
#  if(any(time != as.integer(time))) stop("time must be an integer.")
#  if(any(time<1)) stop("time must be greater than 0.")
#  N = max(time)
#  if(is.null(hyp)){
#    ##Default hyp
#    if(DIST=="Binomial") hyp=c(1,1)	#a,b
#    if(DIST=="Poisson") hyp=c(1,1)	#a,b
#    if(DIST=="Normal") hyp=c(0,1,1,1)	#m,c,a,b
#    if(DIST=="Mean") hyp=c(1,1,1)	#sig2,m,c
#    if(DIST=="Var") hyp=c(0,1,1)	#mu,a,b
#  }
#  if(!is.numeric(hyp)) stop("Invalid prior parameters.")
#  if(DIST=="Poisson" & any(data!=as.integer(data)) & any(data<0)) stop("Invalid data")
#  if(DIST=="Poisson" & length(hyp)!=2 & !all(hyp>0)) stop("Invalid Prior Parameters.")
#  if(DIST=="Binomial" & any(data!=as.integer(data)) & any(data<0) & any(data>1)) stop("Invalid data")
#  if(DIST=="Binomial" & length(hyp)!=2 & !all(hyp>0)) stop("Invalid Prior Parameters.")
#  if(DIST=="Nornal"   & length(hyp)!=4 & !all(hyp[-1]>0)) stop("Invalid Prior Parameters.")
#  if(DIST=="Mean"     & length(hyp)!=3 & !all(hyp>0)) stop("Invalid Prior Parameters.")
#  if(DIST=="Var"      & length(hyp)!=3 & !all(hyp[-1]>0)) stop("Invalid Prior Parameters.")
#
#
#
#  if(!is.numeric(l) | length(l)!=1)  stop("Invalid minimum segment length.")
#  if(l!=as.integer(l) | l<1 | l>N) stop("Invalid minimum segment length.")
#  if(!is.numeric(probs) | any(probs<0) | any(probs>1) | anyNA(probs)) stop("Invalid probs.")
#
#  if(length(cpts)==1) cpts = 1
#  if(is.character(cpts)){
#    HHMM <- seq.POSIXt(from =as.POSIXct("1970-01-01 00:00:00"),
#                       to = as.POSIXct("1970-01-02 00:00:00"),
#                       length.out = N+1)######
#    HHMM <- format(HHMM[-length(HHMM)],"%H:%M")
#    cpts <- sort(which(HHMM %in% cpts))
#  }else{
#    cpts <- sort(cpts)
#  }
#  if(any(cpts<min(time)) | any(cpts>max(time)) | any(diff(c(cpts,cpts[1]+N))<l) )
#    stop("cpts are not a valid vector of changepoint positions.")
#
#  if(length(cpts)<=1){
#    SUM = sum(data)
#    COUNT = length(data)
#    SQ = sum(data^2)
#  }else{
#    SQ <- SUM <- COUNT <- rep(0,length(cpts))
#    L <- time<=cpts[1] | time > cpts[length(cpts)]
#    SUM[1] <- sum(data[L])
#    SQ[1] <- sum(data[L]^2)
#    COUNT[1] <- sum(L)
#    for(j in 2:length(cpts)){
#      L <- time<=cpts[j] & time > cpts[j-1]
#      SUM[j] <- sum(data[L])
#      SQ[j] <- sum(data[L]^2)
#      COUNT[j] <- sum(L)
#    }
#  }
#  out <- matrix(NA,nrow=length(probs),ncol=length(SUM))
#  if(DIST == "Poisson"){
#    ALPHA = hyp[1] + SUM
#    BETA = hyp[2] + COUNT
#    for(p in 1:length(probs)) out[p,] <- qgamma(probs[p],ALPHA,BETA)
#    colnames(out) <- paste0("lam",1:length(SUM))
#  }else if(DIST=="Binomial"){
#    ALPHA = hyp[1] + SUM
#    BETA = hyp[2] + COUNT - SUM
#    for(p in 1:length(probs)) out[p,] <- qbeta(probs[p],ALPHA,BETA)
#    colnames(out) <- paste0("phi",1:length(SUM))
#  }else if(DIST == "Normal"){
#    mean = SUM/COUNT
#    M = (SUM+hyp[1]/hyp[2])/(1/hyp[2]+COUNT)
#    C = 1/(1/hyp[2] + COUNT)
#    A = hyp[3] + COUNT/2
#    B = hyp[4] + 0.5*(SQ - 2*mean*SUM + mean^2) + 0.5*(COUNT/hyp[2])*(mean-hyp[1])^2/(COUNT+1/hyp[2])
#    out2 <- out
#
#    for(p in 1:length(probs)){
#      out2[p,] <- 1/qgamma(1-probs[p],shape=A,rate=B)
#      out[p,] <- qt(probs[p],df=2*A)*sqrt(B*C/A) + M
#    }
#
#    if(length(p)>1){
#      colnames(out) <- paste0("mean",1:length(SUM))
#      colnames(out2) <- paste0("var",1:length(SUM))
#      out <- cbind(out,out2)
#    }else{
#      out <- rbind(out,out2)
#      rownames(out) <- c("mean","var")
#      colnames(out) <- paste0("seg",1:length(SUM))
#    }
#  }else if(DIST == "Mean"){
#    M = (SUM/hyp[1] + hyp[2]/hyp[3])/(COUNT/hyp[1] + 1/hyp[3])
#    C = 1/(1/hyp[3] + COUNT/hyp[1])
#    for(p in 1:length(probs)) out[p,] <- qnorm(probs[p],mean=M,sd=sqrt(C))
#
#  }else if(DIST == "Var"){
#    A = hyp[2] + COUNT/2
#    B = hyp[3] + 0.5*(SQ - 2*hyp[1]*SUM + COUNT*hyp[1]^2)
#    for(p in 1:length(probs)) out[p,] <- 1/qgamma(1-probs[p],shape=A,rate=B)
#  }else{
#    stop("Distribution not recognised.")
#  }
#  if(length(probs)>1){
#    rownames(out) <- paste0(probs*100,"%")
#    return(out)
#  }else{
#    return(out[1,])
#  }
#}



