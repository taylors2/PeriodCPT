setClass("pcpt", slots=list(
    data.set     = "ts",          ## Data set contaning time information
    periodlength = "numeric",     ## Period length (N)
    minseglen    = "numeric",     ## Minimum segment length (L)
    npcpts.max   = "numeric",     ## Maximum number of with period changepoints [floor(N/L)]
    distribution = "character",   ## Distribution
    pcpt.prior   = "list",        ## Prior info for the pcpt prior
    param.prior  = "numeric",     ## Prior info for the segment parameters
    MCMC.options = "list",        ## List of MCMC options
    MCMC.inits   = "list",        ## List of chain initial values
    MCMC.chains  = "list",        ## MCMC output from chain
    param.est    = "list",        ## Mode estimates for segment parameters
    pcpt.est     = "list",        ## Mode estimates for pcpt paramters
    date         = "character",   ## Date that the object was created
    version      = "character"),   ## Version of code

  prototype = prototype(
    date         = date(),
    version      = as(packageVersion("PeriodCPT"), 'character') ) )

############################################################################################
#exportMethods(
#  data.set, periodlength, minseglen, npcpts.max, distribution, pcpt.prior, param.prior,
#  MCMC.options, n.chains, n.iter, n.burn, toggle.quiet, MCMC.inits, MCMC.chain, MCMC.chains,
#  "data.set<-", "periodlength<-", "minseglen<-", "distribution<-", "pcpt.prior<-",
#  "param.prior<-", "MCMC.options<-", "n.chains<-", "n.iter<-", "n.burn<-", "MCMC.inits<-",
#  "MCMC.chain<-", "MCMC.chains<-"
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
setMethod("periodlength","pcpt",function(object) object@periodlength)
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
  if(length(minseglen) != 1 || !is.numeric(minseglen))
    stop("Minimum segment length specifed incorrectly.")
  if(minseglen <= 0 || floor(minseglen) != minseglen)
    stop("Minimum segment length specifed incorrectly.")
  object@minseglen <- minseglen
  if(length(object@periodlength) == 1){
    object@npcpts.max <- floor(object@periodlength / object@minseglen)
  }
  return(minseglen)
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
  if(!any(names(object@MCMC.options) == "n.chains"))
    stop("MCMC.options slot not initialised correctly for `n.chains`.")
  object@MCMC.options$n.chains <- value
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
  if(!any(names(object@MCMC.options) == "n.iter"))
    stop("MCMC.options slot not initialised correctly for `n.iter`.")
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
  if(!any(names(object@MCMC.options) == "n.burn"))
    stop("MCMC.options slot not initialised correctly for `n.burn`.")
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
  if(!any(names(object@MCMC.options) == "quiet"))
    stop("MCMC.options slot not initialised correctly for `quiet`.")
  if(!is.logical(object@MCMC.options$quiet))
    stop("`quiet` item in MCMC.options slot is not logical.")
  object@MCMC.options$quiet <- !object@MCMC.options$quiet})

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

# set/get and reset MCMC slot for given chain index
if(!isGeneric("MCMC.chain")) {
  if(is.function("MCMC.chain")){
    fun <- MCMC.chain
  }else{
    fun <- function(object, index){ standardGeneric("MCMC.chain") }
  }
  setGeneric("MCMC.chain", fun)
}
setMethod("MCMC.chain","pcpt",function(object, index) object@MCMC.chains[[paste0("`",index,"`")]])
setGeneric("MCMC.chain<-", function(object, index, value) standardGeneric("MCMC.chain<-"))
setReplaceMethod("MCMC.chain", "pcpt", function(object, index, value) {
  if(!(index %in% as.numeric(names(object@MCMC))))
    stop(paste0("Index `",index,"` not found in list of chains."))
  object@MCMC.chains[[paste0("`",index,"`")]] <- value
  return(object)
})

# set/get and reset all of MCMC slot
if(!isGeneric("MCMC.chains")) {
  if(is.function("MCMC.chains")){
    fun <- MCMC.chains
  }else{
    fun <- function(object){ standardGeneric("MCMC.chains") }
  }
  setGeneric("MCMC.chains", fun)
}
setMethod("MCMC.chains","pcpt",function(object) object@MCMC.chains)
setGeneric("MCMC.chains<-", function(object, value) standardGeneric("MCMC.chains<-"))
setReplaceMethod("MCMC.chains", "pcpt", function(object, value) {
  object@MCMC.chains <- value
  return(object)
})







