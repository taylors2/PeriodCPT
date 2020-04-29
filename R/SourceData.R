#########################
## GET AND FORMAT DATA ##
#########################
#suppressMessages(library(RPostgreSQL))

hsql <- function(envr,query,config.folder=NULL){
  if(!any(envr == c("live","stage"))) stop("'envr' argument not recognised!")
  #config.folder <-'H:/My Documents/06-Research/R03-Intelesant/'
  if(is.null(config.folder)) stop("Missing config.folder")
  if(max(gregexpr("/",config.folder)[[1]])!=nchar(config.folder)){
    config.folder <- paste0(config.folder,"/")
  }
  config.file <- paste0("config_",envr,".txt")
  file <- paste(config.folder,config.file,sep="")
  if(!file.exists(file)) stop(paste0("File does not exist - ",file))
  config <- read.table(file = file,
                       sep=",",
                       col.names=c("drv", "host", "user", "password", "dbname"),
                       stringsAsFactors=FALSE)
  b <- NULL
######################################
## Need RPostgreSQL package to work ##
######################################
#  con<-dbConnect(drv=config$drv,
#                 host=config$host,
#                 user=config$user,
#                 password=config$password,
#                 dbname=config$dbname)
#  b<-dbGetQuery(conn = con, statement = gsub(pattern='\\n',replacement="",x=query))
#  dbDisconnect(con)
######################################

  return(b)
}

getSensorData <- function(NAME, num.training.days=7, width=15, start=NULL, end=NULL,
                          RAW=FALSE, GroupSensor=TRUE, envr="live",config.folder=getwd()){

  stop("!! This function has stop working from v1.1.5 !!")
  #To re-activate, remove comment marks in 'hsql()' above and 
  #  add 'RPostgreSQL' to package list in Description

  if(length(NAME)>1 || !is.character(NAME)) stop("Invalid NAME")
  if(length(RAW)!=1 || !is.logical(RAW)) stop("Invalid RAW")
  if(length(GroupSensor)!=1 || !is.logical(GroupSensor)) stop("Invalid GroupSensor")
  if(24*60 %% width != 0) stop("width does not factor into minutes in the day.")
  
  #Get time zone information
  TZ <- hsql(envr = envr, config.folder = config.folder,
             query=paste0("SELECT timezone FROM howzsites WHERE hid = '",NAME,"'"))$timezone
  #Format date information
  if(length(start)>1 | length(end)>1) stop("Invalid start and/or end.")
  if(is.null(start)!=is.null(end)) stop("Both start and end must be specified")
  if(!is.null(start)){
    if(!inherits(start,c("character","Date","POSIXct")))
      stop("Invalid start")
    if(!inherits(end,c("character","Date","POSIXct")))
      stop("Invalid end")
    start_date <- as.POSIXct(start, tz=TZ)
    offset <- (3600*as.numeric(format(start_date,"%H")) 
               + 60*as.numeric(format(start_date,"%M")) + as.numeric(format(start_date,"%OS6"))) %% (60*width)
    start_date <- start_date - offset #shift start to beginning of time bin
    start <- as.numeric(start_date)
    end_date <- as.POSIXct(end,tz=TZ)
    offset <- (3600*as.numeric(format(end_date,"%H")) 
               + 60*as.numeric(format(end_date,"%M")) + as.numeric(format(end_date,"%OS6"))) %% (60*width)
    end_date <- end_date - offset - 1 #shift end to last second of previous time bin
    end <- as.numeric(end_date)
    if(start>=end) stop("Start must be earlier than end")
  }else{
    if(length(num.training.days)!=1 || !is.numeric(num.training.days)) stop("Invalid num.training.days")
    if(num.training.days%%1!=0 || num.training.days<=0) stop("Invalid num.training.days")
    
    lasttime <- hsql(envr = envr, config.folder = config.folder,
                     query=paste0("SELECT datecreated FROM messages, messagetext ",
                                  "WHERE messages.textid = messagetext.textid ",
                                  "AND messagetext.activity_tracking = TRUE ",
                                  "AND messages.hid = '",NAME,"' ",
                                  "ORDER BY datecreated DESC "))
    if(nrow(lasttime)==0){
      warning(paste0("Cannot find required data for: ",NAME))
      return(NULL)
    }
    lasttime$time <- as.POSIXct(lasttime$datecreated,origin='1970-01-01',tz=TZ)
    lasttime$date <- as.Date(lasttime$time,tz=TZ)
    tmp <- sort(unique(lasttime$date),decreasing = TRUE)
    
    lasttime_HID <- NA
    if(length(tmp)>num.training.days){
      ind <- 2
      search <- TRUE
      while(search & (ind+num.training.days<=length(tmp))){  
        if(tmp[ind]-num.training.days+1 != tmp[ind+num.training.days-1]){ 
          ind <- ind + 1 
        }else{
          search <- FALSE
        }
      }
      if(!search){
        lasttime_HID <- max(lasttime$time[lasttime$date==tmp[ind]])
      }
    }
    if(is.na(lasttime_HID)){
      warning(paste0("Cannot find required data for: ",NAME))
      return(NULL)
    }
    
    lastobstime <- as.POSIXct(lasttime_HID,origin='1970-01-01',tz=TZ)
    sec <- as.numeric(format(lastobstime,"%OS6"))
    min <- as.numeric(format(lastobstime,"%M"))
    hr  <- as.numeric(format(lastobstime,"%H"))
    secs_of_day <- hr*3600 + min*60 + sec
    offset <- ifelse(secs_of_day>0, 3600*24 - secs_of_day, 0)
    end_date <- lastobstime + offset - 1 
    start_date <- end_date+1 - num.training.days*3600*24
    if(format(start_date,"%H")=="01"){  ##Manually manage 1hr daylight savings if needed
      start_date <- start_date - 3600
    }else if(format(start_date,"%H")=="23"){
      start_date <- start_date + 3600
    }
    end <- as.numeric(end_date)
    start <- as.numeric(start_date)
  }
  data <- hsql(envr = envr, config.folder = config.folder,
               query = paste0("SELECT datecreated, chan, macsidid FROM messages, messagetext ",
                              "WHERE messages.textid = messagetext.textid ",
                              "AND messagetext.activity_tracking = TRUE ",
                              "AND messages.hid = '",NAME,"' ",
                              "AND datecreated BETWEEN ",start," AND ",end,"ORDER BY datecreated"))
  if(nrow(data)==0){
    stop("No data found for date criteria.")
    return(NULL)
  }
  data$datetime <- as.POSIXct(data$datecreated,origin='1970-01-01',tz=TZ)
  if(RAW) return(data)
  data$sensor <- paste0(data$chan,"_",data$macsidid)
  data$binepoch <- format(data$datetime - as.numeric(data$datetime)%%(60*width),"%H:%M")
  
  TBIN <- format(seq.POSIXt(from=start_date, to=end_date, by=60*width),"%Y-%m-%d %H:%M")
  
  data$binepoch <- paste0(as.character(as.Date(data$datetime, tz=TZ))," ",data$binepoch)
  SensorBIN <- unique(data$sensor)
  BINCOUNT <- BINACTIVE <- matrix(0,nrow=length(TBIN),ncol=length(SensorBIN))
  for(s in 1:length(SensorBIN)){
    for(i in 1:length(TBIN)){
      L <- data$sensor==SensorBIN[s] & data$binepoch==TBIN[i]
      BINCOUNT[i,s] <- sum(L)
      BINACTIVE[i,s] <- as.numeric(any(L))
    }
  }
  
  if(GroupSensor){
    BinCount <- apply(BINCOUNT,1,sum)
    BinActive <- as.numeric(apply(BINACTIVE==1,1,any))
    out <- data.frame(Time=TBIN,Count=BinCount,Active=BinActive)
  }else{
    BinCount <- as.data.frame(BINCOUNT)
    names(BinCount) <- paste0("Count_",SensorBIN)
    BinActive <- as.data.frame(BINACTIVE)
    names(BinActive) <- paste0("Active_",SensorBIN)
    out <- cbind(data.frame(Time=TBIN),BinCount,BinActive)
  }
  out$binid <- unlist(lapply(lapply(strsplit(as.character(out$Time)," "),
                                    function(a){a[2]}),function(b){sum(as.numeric(unlist(
                                      strsplit(b,":")))*c(60,1))}))/width + 1
  return(out)
}
