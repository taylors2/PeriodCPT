setClass("pcpt", slots=list(
    data.set     = "ts",          ## Data set contaning time information
    periodlength = "numeric",     ## Period length (N)
    minseglen    = "numeric",     ## Minimum segment length (L)
    npcpts.max   = "numeric",     ## Maximum number of with period changepoints [floor(N/L)]
    distribution = "character",   ## Distribution
    pcpt.prior   = "list",        ## Prior info for the pcpt prior
    param.prior  = "numeric",     ## Prior info for the segment parameters
    MCMC.options = "list",        ## List of MCMC options
    MCMC         = "list",        ## MCMC output from chain
    param.est    = "list",        ## Mode estimates for segment parameters
    pcpt.est     = "list",        ## Mode estimates for pcpt paramters
    date         = "character",   ## Date that the object was created
    version      = "character")  ## Version of code
)
#  prototype = prototype(
#    date         = date(),
#    version      = as(packageVersion("PeriodCPT"), 'character') ) )

############################################################################################
#data.set       data.set<-
#periodlength   periodlength<-
#minseglen      minseglen<-
#npcpts.max
#distribution   distribution<-
#pcpt.prior     pcpt.prior<-
#param.prior    param.prior<-
#MCMC.options   MCMC.options<-
  #n.chains
  #n.iter         n.iter<-
  #n.burn         n.burn<-
  #toggle.quiet


{

## set/get and reset data.set slot
if(!isGeneric("data.set")) {
  if (is.function("data.set")){
    fun <- data.set
  }else {
    fun <- function(object){ standardGeneric("data.set") }
  }
  setGeneric("data.set", fun)
}
setMethod("data.set","pcpt",function(object) coredata(object@data.set))
setGeneric("data.set<-", function(object, value) standardGeneric("data.set<-"))
setReplaceMethod("data.set", "pcpt", function(object, value) {
  ##Input check
  if(!is.numeric(value)) stop("Data must contain numeric values.")
  if(anyNA(value)) stop("Data must not contain missing values.")
  if(is.ts(value)){
    object@data.set <- value
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
setMethod("periodlength","pcpt",function(object) coredata(object@periodlength))
setGeneric("periodlength<-", function(object, value) standardGeneric("periodlength<-"))
setReplaceMethod("periodlength", "pcpt", function(object, value) {
  if(is.null(periodlength)){
    object@periodlength <- round(1/frequency(data.set(object)))
  }else{
    if(length(periodlength) != 1 || !is.numeric(periodlength))
      stop("Period length specified incorrectly.")
    if(floor(periodlength) != periodlength || periodlength < 0)
      stop("Period length specified incorrectly.")
    object@periodlength <- periodlength
  }
  if(length(object@minseglen) == 1){
    object@npcpts.max <- floor(object@periodlength / object@minseglen)
  }
  return(object)
}

## set/get and reset minseglen slot
if(!isGeneric("minseglen")) {
  if (is.function("minseglen")){
    fun <- minseglen
  }else {
    fun <- function(object){ standardGeneric("minseglen") }
  }
  setGeneric("minseglen", fun)
}
setMethod("minseglen","pcpt",function(object) coredata(object@minseglen))
setGeneric("minseglen<-", function(object, value) standardGeneric("minseglen<-"))
setReplaceMethod("minseglen", "pcpt", function(object, value) {
  if(length(minseglen) != 1 || !is.numeric(minseglen))
    stop("Minimum segment length specifed incorrectly.")
  if(minseglen < 0 || floor(minseglen) != minseglen)
    stop("Minimum segment length specifed incorrectly.")
  object@minseglen <- minseglen
  if(length(object@periodlength) == 1){
    object@npcpts.max <- floor(object@periodlength / object@minseglen)
  }
  return(minseglen)
}

## set/get npcpts.max slot
if(!isGeneric("npcpts.max")) {
  if (is.function("npcpts.max")){
    fun <- npcpts.max
  }else {
    fun <- function(object){ standardGeneric("npcpts.max") }
  }
  setGeneric("npcpts.max", fun)
}
setMethod("npcpts.max","pcpt",function(object) coredata(object@npcpts.max))

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
  object@MCMC.options$n.iter <- value
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
  object@MCMC.options$n.burn <- value
  return(object)
})

if(!isGeneric("toggle.quiet")) {
  if(is.function("toggle.quiet")){
    fun <- toggle.quiet
  }else{
    fun <- function(object){ standardGeneric("toggle.quiet") }
  }
  setGeneric("toggle.quiet", fun)
}
setMethod("toggle.quiet","pcpt",function(object){
  object@MCMC.options$quiet <- !object@MCMC.options$quiet})

}




###############################################################
#tmp <- new("pcpt")
#data.set(tmp)
#tmpdata <- ts(rnorm(10))
#data.set(tmp) <- tmpdata
#data.set(tmp)
#tmp



