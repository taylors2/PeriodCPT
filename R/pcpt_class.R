setClass("pcpt", slots=list(
    data.set     = "ts",          ## Data set contaning time information
    periodlength = "numeric",     ## Period length (N)
    minseglen    = "numeric",     ## Minimum segment length (L)
    npcpts.max   = "numeric",     ## Maximum number of with period changepoints [floor(N/L)]
    distribution = "character",   ## Distribution
    nsegparam    = "numeric",     ## Number of paramters per segment.
    pcpt.prior   = "list",        ## Prior info for the pcpt prior
    param.prior  = "numeric",     ## Prior info for the segment parameters
    MCMC.options = "list",        ## List of MCMC options
    MCMC.inits   = "list",        ## List of chain initial values
    MCMC.last    = "list",        ## List of chain values at last interation (filled at summarised)
    results      = "list",        ## MCMC output from chain / Summarised chain(s)
    summarised   = "logical",     ## Does results slot are the raw chain (FALSE) or are summaised (TRUE)
    param.mode   = "matrix",      ## Mode estimates for segment parameters
    pcpt.mode    = "numeric",     ## Mode estimates for pcpt paramters
    fit          = "numeric",     ## Fit metric, evaluated at summary
    date         = "character",   ## Date that the object was created
    version      = "character"),  ## Version of code

  prototype = prototype(
    version      = as(packageVersion("PeriodCPT"), 'character') ) )

############################################################################################
#exportMethods(
#  data.set, periodlength, minseglen, npcpts.max, distribution, pcpt.prior, param.prior,
#  MCMC.options, n.chains, n.iter, n.burn, toggle.quiet, MCMC.inits, result, results,
#  "data.set<-", "periodlength<-", "minseglen<-", "distribution<-", "pcpt.prior<-",
#  "param.prior<-", "MCMC.options<-", "n.chains<-", "n.iter<-", "n.burn<-", "MCMC.inits<-",
#  "result<-", "results<-"
#)

##TODO!!!
#########
#nseg/ncpts
#cpts
#param.est/param/coef
#seg.len
#plot
#summary
#show
#print
#logLik/likelihood
#quantile
############################################################################################



## set/get and reset data.set slot
if(!isGeneric("data.set")) {
  if (is.function("data.set")){
    fun <- data.set
  }else {
    fun <- function(object){ standardGeneric("data.set") }
  }
  setGeneric("data.set", fun)
}
setMethod("data.set","pcpt",function(object) object@data.set)
setGeneric("data.set<-", function(object, value) standardGeneric("data.set<-"))
setReplaceMethod("data.set", "pcpt", function(object, value) {
  ##Input check
  if(is.null(value)) stop("Data is missing.")
  if(anyNA(value)) stop("Data must not contain NA missing values.")
  if(is.ts(value)){
    object@data.set <- value
    periodlength(object) <- frequency(value)
  }else{
    object@data.set <- ts(value)
  }
  return(object)
})

## set/get and reset periodlength slot
if(!isGeneric("periodlength")) {
  if (is.function("periodlength")){
    fun <- periodlength
  }else {
    fun <- function(object){ standardGeneric("periodlength") }
  }
  setGeneric("periodlength", fun)
}
setMethod("periodlength","pcpt",function(object) object@periodlength)
setGeneric("periodlength<-", function(object, value) standardGeneric("periodlength<-"))
setReplaceMethod("periodlength", "pcpt", function(object, value) {
  if(is.null(value) & is.ts(data.set(object))){
    value <- frequency(data.set(object))
  }
  if(length(value) != 1 || !is.numeric(value))
    stop("Period length specified incorrectly.")
  if(floor(value) != value || value <= 0)
    stop("Period length specified incorrectly.")
  #if(value == 1) stop("Period length must be greater than 1.")
  object@periodlength <- value

  if(length(object@minseglen) == 1){
    if(object@minseglen > object@periodlength){
      stop("Minimum segment length longer than period length.")
    }
    object@npcpts.max <- floor(object@periodlength / object@minseglen)
  }
  return(object)
})

## set/get and reset minseglen slot
if(!isGeneric("minseglen")) {
  if (is.function("minseglen")){
    fun <- minseglen
  }else {
    fun <- function(object){ standardGeneric("minseglen") }
  }
  setGeneric("minseglen", fun)
}
setMethod("minseglen","pcpt",function(object) object@minseglen)
setGeneric("minseglen<-", function(object, value) standardGeneric("minseglen<-"))
setReplaceMethod("minseglen", "pcpt", function(object, value) {
  if(length(value) != 1 || !is.numeric(value))
    stop("Minimum segment length specifed incorrectly.")
  if(value <= 0 || floor(value) != value)
    stop("Minimum segment length specifed incorrectly.")
  object@minseglen <- value
  if(length(object@periodlength) == 1){
    if(value > periodlength(object)){
      stop("Minimum segment length longer than period length.")
    }
    object@npcpts.max <- floor(object@periodlength / object@minseglen)
  }
  return(object)
})

## set/get npcpts.max slot
if(!isGeneric("npcpts.max")) {
  if (is.function("npcpts.max")){
    fun <- npcpts.max
  }else {
    fun <- function(object){ standardGeneric("npcpts.max") }
  }
  setGeneric("npcpts.max", fun)
}
setMethod("npcpts.max","pcpt",function(object) object@npcpts.max)

# set/get and reset distribution slot
if(!isGeneric("distribution")) {
  if(is.function("distribution")){
    fun <- distribution
  }else{
    fun <- function(object){ standardGeneric("distribution") }
  }
  setGeneric("distribution", fun)
}
setMethod("distribution","pcpt",function(object) object@distribution)
setGeneric("distribution<-", function(object, value) standardGeneric("distribution<-"))
setReplaceMethod("distribution", "pcpt", function(object, value) {
  object@distribution <- value
  return(object)
})

# set/get and reset pcpt.prior slot
if(!isGeneric("pcpt.prior")) {
  if(is.function("pcpt.prior")){
    fun <- pcpt.prior
  }else{
    fun <- function(object){ standardGeneric("pcpt.prior") }
  }
  setGeneric("pcpt.prior", fun)
}
setMethod("pcpt.prior","pcpt",function(object) object@pcpt.prior)
setGeneric("pcpt.prior<-", function(object, value) standardGeneric("pcpt.prior<-"))
setReplaceMethod("pcpt.prior", "pcpt", function(object, value) {
  object@pcpt.prior <- value
  return(object)
})

# set/get and reset param.prior slot
if(!isGeneric("param.prior")) {
  if(is.function("param.prior")){
    fun <- param.prior
  }else{
    fun <- function(object){ standardGeneric("param.prior") }
  }
  setGeneric("param.prior", fun)
}
setMethod("param.prior","pcpt",function(object) object@param.prior)
setGeneric("param.prior<-", function(object, value) standardGeneric("param.prior<-"))
setReplaceMethod("param.prior", "pcpt", function(object, value) {
  object@param.prior <- value
  return(object)
})


# set/get and reset MCMC.options slot
if(!isGeneric("MCMC.options")) {
  if(is.function("MCMC.options")){
    fun <- MCMC.options
  }else{
    fun <- function(object){ standardGeneric("MCMC.options") }
  }
  setGeneric("MCMC.options", fun)
}
setMethod("MCMC.options","pcpt",function(object) object@MCMC.options)
setGeneric("MCMC.options<-", function(object, value) standardGeneric("MCMC.options<-"))
setReplaceMethod("MCMC.options", "pcpt", function(object, value) {
  object@MCMC.options <- value
  return(object)
})

# set/get n.chains item in MCMC.options slot
if(!isGeneric("n.chains")) {
  if(is.function("n.chains")){
    fun <- n.chains
  }else{
    fun <- function(object){ standardGeneric("n.chains") }
  }
  setGeneric("n.chains", fun)
}
setMethod("n.chains","pcpt",function(object) object@MCMC.options$n.chains)
setGeneric("n.chains<-", function(object, value) standardGeneric("n.chains<-"))
setReplaceMethod("n.chains", "pcpt", function(object, value) {
  object@MCMC.options[["n.chains"]] <- value
  return(object)
})

# set/get and reset n.iter item in MCMC.options slot
if(!isGeneric("n.iter")) {
  if(is.function("n.iter")){
    fun <- n.iter
  }else{
    fun <- function(object){ standardGeneric("n.iter") }
  }
  setGeneric("n.iter", fun)
}
setMethod("n.iter","pcpt",function(object) object@MCMC.options$n.iter)
setGeneric("n.iter<-", function(object, value) standardGeneric("n.iter<-"))
setReplaceMethod("n.iter", "pcpt", function(object, value) {
  object@MCMC.options[["n.iter"]] <- value
  return(object)
})

# set/get and reset n.burn item in MCMC.options slot
if(!isGeneric("n.burn")) {
  if(is.function("n.burn")){
    fun <- n.burn
  }else{
    fun <- function(object){ standardGeneric("n.burn") }
  }
  setGeneric("n.burn", fun)
}
setMethod("n.burn","pcpt",function(object) object@MCMC.options$n.burn)
setGeneric("n.burn<-", function(object, value) standardGeneric("n.burn<-"))
setReplaceMethod("n.burn", "pcpt", function(object, value) {
  object@MCMC.options[["n.burn"]] <- value
  return(object)
})

if(!isGeneric("quiet")) {
  if(is.function("quiet")){
    fun <- quiet
  }else{
    fun <- function(object){ standardGeneric("quiet") }
  }
  setGeneric("quiet", fun)
}
setMethod("quiet","pcpt",function(object) object@MCMC.options$quiet)
setGeneric("quiet<-", function(object, value) standardGeneric("quiet<-"))
setReplaceMethod("quiet", "pcpt", function(object, value) {
  object@MCMC.options[["quiet"]] <- value
  return(object)
})


# set/get and reset MCMC.inits slot
if(!isGeneric("MCMC.inits")) {
  if(is.function("MCMC.inits")){
    fun <- MCMC.inits
  }else{
    fun <- function(object){ standardGeneric("MCMC.inits") }
  }
  setGeneric("MCMC.inits", fun)
}
setMethod("MCMC.inits","pcpt",function(object) object@MCMC.inits)
setGeneric("MCMC.inits<-", function(object, value) standardGeneric("MCMC.inits<-"))
setReplaceMethod("MCMC.inits", "pcpt", function(object, value) {
  object@MCMC.inits <- value
  return(object)
})


# set/get and reset MCMC.last slot
if(!isGeneric("MCMC.last")) {
  if(is.function("MCMC.last")){
    fun <- MCMC.last
  }else{
    fun <- function(object){ standardGeneric("MCMC.last") }
  }
  setGeneric("MCMC.last", fun)
}
setMethod("MCMC.last","pcpt",function(object) object@MCMC.last)
setGeneric("MCMC.last<-", function(object, value) standardGeneric("MCMC.last<-"))
setReplaceMethod("MCMC.last", "pcpt", function(object, value) {
  object@MCMC.last <- value
  return(object)
})

# set/get and reset MCMC slot for given chain index
if(!isGeneric("result")) {
  if(is.function("result")){
    fun <- result
  }else{
    fun <- function(object, index){ standardGeneric("result") }
  }
  setGeneric("result", fun)
}
setMethod("result","pcpt",function(object, index){
  if(missing(index)) stop("Argument 'index' missing with no default.")
  if(!is.numeric(index) | length(index) != 1) stop("Argument 'index' must be a single numeric value.")
  if(!(index %in% as.numeric(names(object@results)))){
    if((index > length(summarised(object))) & (index <= n.chains(object))){
      stop(paste0("Cannot access information for index `",index,
                  "`. It is likely that the contents accoss multiple chains ",
                  "has been combined into the first index."))
    }else{
      stop(paste0("Index `",index,"` not found in results list."))
    }
  }
  object@results[[as.character(index)]]
})
setGeneric("result<-", function(object, index, value) standardGeneric("result<-"))
setReplaceMethod("result", "pcpt", function(object, index, value) {
  if(!(index %in% as.numeric(names(object@results))))
    stop(paste0("Index `",index,"` not found in results list."))
  object@results[[as.character(index)]] <- value
  return(object)
})

# set/get and reset all of MCMC slot
if(!isGeneric("results")) {
  if(is.function("results")){
    fun <- results
  }else{
    fun <- function(object){ standardGeneric("results") }
  }
  setGeneric("results", fun)
}
setMethod("results","pcpt",function(object) object@results)
setGeneric("results<-", function(object, value) standardGeneric("results<-"))
setReplaceMethod("results", "pcpt", function(object, value) {
  object@results <- value
  return(object)
})

if(!isGeneric("summarised")) {
  if(is.function("summarised")){
    fun <- summarised
  }else{
    fun <- function(object){ standardGeneric("summarised") }
  }
  setGeneric("summarised", fun)
}
setMethod("summarised","pcpt",function(object) object@summarised)
setGeneric("summarised<-", function(object, value) standardGeneric("summarised<-"))
setReplaceMethod("summarised", "pcpt", function(object, value) {
  if(!is.logical(value) | anyNA(value))
    stop("Can only assign logical to summarised slot.")
  object@summarised <- value
  return(object)
})

if(!isGeneric("summarized")) {
  if(is.function("summarized")){
    fun <- summarized
  }else{
    fun <- function(object){ standardGeneric("summarized") }
  }
  setGeneric("summarized", fun)
}
setMethod("summarized","pcpt",function(object) summarised(object))
setGeneric("summarized<-", function(object, value) standardGeneric("summarized<-"))
setReplaceMethod("summarized", "pcpt", function(object, value) {
  summarised(object) <- value
  return(object)
})

if(!isGeneric("nsegparam")) {
  if(is.function("nsegparam")){
    fun <- nsegparam
  }else{
    fun <- function(object){ standardGeneric("nsegparam") }
  }
  setGeneric("nsegparam", fun)
}
setMethod("nsegparam","pcpt",function(object) object@nsegparam)
setGeneric("nsegparam<-", function(object, value) standardGeneric("nsegparam<-"))
setReplaceMethod("nsegparam", "pcpt", function(object, value) {
  if(length(value)!=1 | !is.numeric(value) | any(value<1) | anyNA(value)){
    stop("Assignment to nsegparam slot is not a single positive integer.")
  }else if(any(value != floor(value))){
    stop("Assignment to nsegparam slot is not a single positive integer.")
  }
  object@nsegparam <- value
  return(object)
})



#################################################

setMethod("show","pcpt",function(object){
  cat("Class 'pcpt' : Changepoint Object\n")
  cat("        ~~   : S4 class containing", length(attributes(object))-1, "slots with names\n")
  cat("             ", names(attributes(object))[1:(length(attributes(object))-1)], "\n\n")
  cat("Created on   :", object@date, "\n\n")
  cat("summary(.)   :\n----------\n")
  summary(object)
})

setMethod("print","pcpt",function(x, ...){
  show(x)
})

setMethod("summary","pcpt",function(object, ...){

  pcpt_mode  <- pcpt.mode(object)
  param_mode <- param.mode(object)
  cat("Created Using changepoint version",object@version,'\n')
  cat("Distribution            : ", distribution(object), '\n')
  cat("Period length           : ", periodlength(object), "\n")
  cat("Minimum Segment Length  : ", minseglen(object),    "\n")
  cat("Maximum no. of cpts     : ", npcpts.max(object),   "\n")
  cat("Number of chains        : ", n.chains(object),     "\n")
  cat("Number of periodic segs : ", nsegs(object),        "\n")
  if(anyNA(nsegs(object))){
    cat("Periodic cpt locations  : NA\n")
  }else if(nsegs(object) == 1){
    cat("Periodic cpt locations  : (null)\n")
  }else{
    tau <- paste0(unname(as.numeric(pcpt.mode(object))), collapse = ", ")
    cat("Periodic cpt locations  : ", tau,"\n")
  }
  if(anyNA(nsegs(object))){
    cat("Seg. parameters at mode : NA\n")
  }else{
    cat("Seg. parameters at mode : \n")
    print(param.mode(object))
  }
})

setMethod("quantile","pcpt",function(x, ...){
  return(quantile.master(object = x, ...))
})

setMethod("plot","pcpt",function(x, ...){
  cat("Plot function for pcpt class not yet implemented !!!\n")
})


##################################################################
# set/get and reset pcpt.mode slot
if(!isGeneric("pcpt.mode")) {
  if(is.function("pcpt.mode")){
    fun <- pcpt.mode
  }else{
    fun <- function(object){ standardGeneric("pcpt.mode") }
  }
  setGeneric("pcpt.mode", fun)
}
setMethod("pcpt.mode","pcpt",function(object) object@pcpt.mode)
setGeneric("pcpt.mode<-", function(object, value) standardGeneric("pcpt.mode<-"))
setReplaceMethod("pcpt.mode", "pcpt", function(object, value) {
  object@pcpt.mode <- value
  return(object)
})

# set/get and reset param.mode slot
if(!isGeneric("param.mode")) {
  if(is.function("param.mode")){
    fun <- param.mode
  }else{
    fun <- function(object){ standardGeneric("param.mode") }
  }
  setGeneric("param.mode", fun)
}
setMethod("param.mode","pcpt",function(object) object@param.mode)
setGeneric("param.mode<-", function(object, value) standardGeneric("param.mode<-"))
setReplaceMethod("param.mode", "pcpt", function(object, value) {
  object@param.mode <- value
  return(object)
})

# set/get and reset fit slot for given chain index
if(!isGeneric("fit")) {
  if(is.function("fit")){
    fun <- fit
  }else{
    fun <- function(object){ standardGeneric("fit") }
  }
  setGeneric("fit", fun)
}
setMethod("fit","pcpt",function(object) object@fit)
setGeneric("fit<-", function(object, value) standardGeneric("fit<-"))
setReplaceMethod("fit", "pcpt", function(object, value) {
  object@fit <- value
  return(object)
})


## set/get and reset data.set slot
if(!isGeneric("seglen")) {
  if (is.function("seglen")){
    fun <- seglen
  }else {
    fun <- function(object){ standardGeneric("seglen") }
  }
  setGeneric("seglen", fun)
}
setMethod("seglen","pcpt",function(object){
  if(length(pcpt.mode(object))==0){
    warning("Mode pcpt has not yet been evaluated.")
  }else{
    tau <- pcpt.mode(object)
    sl <- diff(c(tau[length(tau)]-periodlength(object), tau))
    unname(sl)
    return(sl)
  }
})


## set/get and reset data.set slot
if(!isGeneric("nsegs")) {
  if (is.function("nsegs")){
    fun <- nsegs
  }else {
    fun <- function(object){ standardGeneric("nsegs") }
  }
  setGeneric("nsegs", fun)
}
setMethod("nsegs","pcpt",function(object){
  if(length(pcpt.mode(object))==0){
    warning("Mode pcpt has not yet been evaluated.")
    return(NA)
  }else{
    return(length(pcpt.mode(object)))
  }
})

## set/get and reset data.set slot
if(!isGeneric("npcpts")) {
  if (is.function("npcpts")){
    fun <- npcpts
  }else {
    fun <- function(object){ standardGeneric("npcpts") }
  }
  setGeneric("npcpts", fun)
}
setMethod("npcpts","pcpt",function(object){return(nsegs(object))})


