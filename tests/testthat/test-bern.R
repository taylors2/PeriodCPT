##Function to generate test data based on selected case
make_test_data <- function(dist, caseid){
  set.seed(1)
  periodlength <- 24
  cycles <- 10

  if(dist == "bern"){
    if(caseid == 1){
      data <- ts(data = rbinom(cycles*periodlength, size = 1,
                               prob = rep( rep( c(0.25, 0.75), each = periodlength/2), cycles) ),
                 frequency = periodlength)
      return(data)
    }
    data  <- ts(data = rbinom(cycles*periodlength, size = 1, prob = 0.25),
                     frequency = periodlength)
    if(caseid == 2) return(data)
  }else if(dist == "pois"){
    if(caseid == 1){
      data <- ts(data = rpois(cycles*periodlength,
                              lambda = rep( rep( c(1, 4), each = periodlength/2), cycles) ),
                              frequency = periodlength)
      return(data)
    }
    data  <- ts(data = rpois(cycles*periodlength, lambda = 1),
                            frequency = periodlength)
    if(caseid == 2) return(data)
  }else if(dist == "mean"){
    if(caseid == 1){
      data <- ts(data = rnorm(cycles*periodlength, sd = 1,
                   mean = rep( rep( c(5, -5), each = periodlength/2), cycles) ),
                   frequency = periodlength)
      return(data)
    }
    data  <- ts(data = rnorm(cycles*periodlength, mean = 5, sd = 1),
                            frequency = periodlength)
    if(caseid == 2) return(data)
  }else if(dist == "var"){
    if(caseid == 1){
      data <- ts(data = rnorm(cycles*periodlength, mean = 0,
                    sd = rep( rep( c(1, 3), each = periodlength/2), cycles) ),
                   frequency = periodlength)
      return(data)
    }
    data <- ts(data = rnorm(cycles*periodlength, mean = 0, sd = 1),
                           frequency = periodlength)
    if(caseid == 2) return(data)
  }else if(dist == "norm"){
    if(caseid == 1){
      data <- ts(data = rnorm(cycles*periodlength,
                   mean = rep( rep( c(-5, 5), each = periodlength/2), cycles),
                   sd = rep( rep( c( 1, 3), each = periodlength/2), cycles) ),
                  frequency = periodlength)
      return(data)
    }
    data <- ts(data = rnorm(cycles*periodlength, mean = 5, sd = 1),
                            frequency = periodlength)
    if(caseid == 2) return(data)
  }else{
    data <- ts(rep(0, cycles*periodlength), frequency = periodlength)
    return(data)
  }

  if(caseid == 3){   ##zeros
    data   <- ts(rep(0, cycles*periodlength), frequency = periodlength)
  }else if(caseid == 4){  ##NAs
    data[sample.int(length(data), size = 10)] <- NA
  }else if(caseid == 5){  ##negative
    data <- -abs(data)
  }else if(caseid == 6){  ##non-integer
    data <- round(data)+0.5
  }else if(caseid == 7){  ##non-binary
    data <- (data < mean(data)) + 1
  }else if(caseid == 8){  ##short
    data <- c(0)
  }else if(caseid == 9){  ##non-ts
    data <- as.numeric(data)
  }else if(caseid == 10){  #char
    data <- ts(c("A","B")[(data>mean(data)) + 1], frequency = periodlength)
  }else if(caseid == 11){  #factor
    data <- ts(as.factor(c("A","B")[(data>mean(data)) + 1]), frequency = periodlength)
  }else if(caseid == 12){ ##non-integer data frequency
    data <- ts(data, frequency = 5.5)
  }else if(caseid == 13){ ##frequency of 1 (ie no period)
    data <- ts(data, frequency = 1)
  }else if(caseid == 14){
    data <- ts(cbind(data,data), frequency = frequency(data))
  }else{  #NULL
    data <- NULL
  }
  return(data)
}

#Lists defining all options per input argument
##############################################

options_distribution <- list("bern","pois","mean","var","norm","invalid")

options_periodlength <- list(24, #valid
                             NULL, #Invalid only if data is not ts
                             1, NA, -1, 24.5, 241, "A") #Invalid
options_minseglen <- list( 2, #valid
                           NA, NULL, -1, 0.1, 1+24, "A") #invalid
options_Mprior <- list("pois","unif",          #valid
                       "invalid", 1, NA, NULL) #invalid
options_Mhyp <- list(1, #valid
                     NA, NULL, 0, -1, "A") #invalid    (nb: only applicable when Mprior == "pois")
options_spread <- list(1, #valid
                       NA, NULL, 0, -1, "A") #invalid
options_param.a <- options_param.b <- options_param.m <-
  options_param.c <- options_const.var <- list(1, #valid
                                           NA, NULL, 0, -1, "A") #invalid
valid_inits_fn1   <- function( ... ) return(1)
valid_inits_fn2   <- function( pcpt.object, chain.index, ... ) return(1)
invalid_inits_fn1 <- function(unspecified_input) return(0)
invalid_inits_fn2 <- function(pcpt.object) return(0)
invalid_inits_fn3 <- function(chain.index) return(0)

options_inits <- list(
  NULL,                            ##should simulate from prior
  "ends",                          ##should over-ride n.chain
  valid_inits_fn1,                 ##should define value from function
  valid_inits_fn2,                 ##should define value from function
  list(1,1),                       #Valid if n.chain == 2
  c(1,2),                          #valid if minseglen == 1
  invalid_inits_fn1,               #invalid
  invalid_inits_fn2,               #invalid
  invalid_inits_fn3,               #invalid
  c(0,12/2),             #invalid
  c(1,12 + 1),           #invalid
  c(0,12/2),             #invalid
  c(12,12/2),  #invalid
  c(1,12),               #invalid
  1:(12+1)               #invalid
)

options_DOTDOTDOT <- list(inits_fn_args = 0,  #valid
                          pcpt.object=new("pcpt"), chain.index=1) #invalid

options_n.chains <- list(1,2, #valid
                        NA,NULL,0,-1,"A",1.1) #invlaid
options_n.iter <- list(100,   #valid
                      NA,NULL,0,-1,"A",1.5) #invalid
options_n.burn <- list(100,0,   #valid
                      NA,NULL,-1,"A",5.5) #invalid
options_cachesize <- list(50,  #valid
                         NA,NULL,0,-1,"A",50.5) #invalid
options_quiet <- list(TRUE, FALSE,  #valid
                     NA,NULL,0,-1,"A",0.5) #invalid


##Error messages:
# statements like $(sub_st)$options_Mprior[[case[5]]]$(sub_ed)$ are used for general messages where an appropriate
#    substitution is made for the given case.
ErrorMessages <- c(
  "gettextf('%s should be one of %s', \"'arg'\", paste(dQuote(c('bern', 'pois', 'norm', 'mean', 'var')), collapse = ', '))",
  "gettextf('%s should be one of %s', \"'arg'\", paste(dQuote(c('pois','unif')), collapse = ', '))",
  "Data is missing.",
  "Data must not contain NA missing values.",
  "Period length is not defined either via data as `ts` object or explicitly given as input.",
  "Mhyp specified incorrectly for Mprior `$(sub_st)$options_Mprior[[case[5]]]$(sub_ed)$`.",
  "Implementation Error: Mprior `$(sub_st)$options_Mprior[[case[5]]]$(sub_ed)$` is not supported.",
  "Hyper-parameter `spread` specified incorrectly.",
  "MCMC option - n.iter specified incorrectly.",
  "MCMC option - n.chains specified incorrectly.",    ###10
  "MCMC option - n.burn specified incorrectly.",
  "MCMC option - cachesize specified incorrectly.",
  "MCMC option - quiet specified incorrectly.",
  "Unrecognised prior for the number of within period changepoints.",
  "Cannot pass arguments `pcpt.object` and `chain.index` from `...` into inits(). These inputs to inits() are managed internally.",
  "Too few inital values provided for specified number of chains.",
  "Incorrect number of initial values for specified number of chains.",
  "Class of at least one inits is not numeric or interger.",
  "Incorrect number of within period changepoints specified by inits.",
  "In inits, within period cpts must be whole numbers.",    ###20
  "In inits, within period cpts must be within [1, period length].",
  "In inits, within period cpts does not satisfy minimum segment length condition.",
  "Hyper-parameter `param.a` specified incorrectly.",
  "Hyper-parameter `param.b` specified incorrectly.",
  "Hyper-parameter `param.c` specified incorrectly.",
  "Hyper-parameter `param.m` specified incorrectly.",
  "Known constant `const.var` is specified incorrectly.",
  "Data is invalid for '$(sub_st)$options_distribution[[case[1]]]$(sub_ed)$' sampling distribution.",
  "Data must contain numeric values.",
  "Cannot find and assign period length.",    ###30
  "Data frequency must be an integer.",
  "Period length specified incorrectly.",
  "Minimum segment length specifed incorrectly.",
  "MCMC.options slot not initialised correctly for `n.chains`.",
  "MCMC.options slot not initialised correctly for `n.iter`.",
  "MCMC.options slot not initialised correctly for `n.burn`.",
  "MCMC.options slot not initialised correctly for `quiet`.",
  "Index `$(sub_st-ish)$index$(sub_ed-ish)$` not found in list.",   ##I have no idea how to invoke this error message!!!
  "Can only assign logical to summarised slot.",
  "Assignment to nseglen slot is not a single positive integer.",    ###40
  "Length of data is too short or period length is too long.",
  "Period length must be greater than 1.",
  "Minimum segment length longer than period length.",
  "'arg' must be NULL or a character vector",
  "Mprior cannot be NULL.",
  "MCMC option - n.iter not specified."
  )

##Function to perform the testthat commands
PeriodCPT_TEST <- function(case){
  ##Generate data
  data <- make_test_data(options_distribution[[case[1]]], case[2])
  if(case[length(case)]==0){
      #Expect function to be valid, returning a "pcpt" class object
      expect_s4_class(PeriodCPT(data = data,
                          distribution = options_distribution[[ case[01] ]],
                          periodlength = options_periodlength[[ case[03] ]],
                          minseglen    = options_minseglen[[    case[04] ]],
                          Mprior       = options_Mprior[[       case[05] ]],
                          Mhyp         = options_Mhyp[[         case[06] ]],
                          spread       = options_spread[[       case[07] ]],
                          param.a      = options_param.a[[      case[08] ]],
                          param.b      = options_param.b[[      case[09] ]],
                          param.c      = options_param.c[[      case[10] ]],
                          param.m      = options_param.m[[      case[11] ]],
                          const.var    = options_const.var[[    case[12] ]],
                          inits        = options_inits[[        case[13] ]],
                          n.iter       = options_n.iter[[       case[14] ]],
                          n.chains     = options_n.chains[[     case[15] ]],
                          n.burn       = options_n.burn[[       case[16] ]],
                          cachesize    = options_cachesize[[    case[17] ]],
                          quiet        = options_quiet[[        case[18] ]]),
                "pcpt")
  }else{
      #Expect test to return error message

      #Format error message for case
      msg <- ErrorMessages[case[length(case)]]
      msg <- gsub('$(sub_st)$', '",', msg, fixed=TRUE)
      msg <- gsub('$(sub_ed)$', ',"', msg, fixed=TRUE)
      msg <- eval(parse(text = paste0('paste0("',msg,'")')))

      expect_that(PeriodCPT(data = data,
                          distribution = options_distribution[[ case[01] ]],
                          periodlength = options_periodlength[[ case[03] ]],
                          minseglen    = options_minseglen[[    case[04] ]],
                          Mprior       = options_Mprior[[       case[05] ]],
                          Mhyp         = options_Mhyp[[         case[06] ]],
                          spread       = options_spread[[       case[07] ]],
                          param.a      = options_param.a[[      case[08] ]],
                          param.b      = options_param.b[[      case[09] ]],
                          param.c      = options_param.c[[      case[10] ]],
                          param.m      = options_param.m[[      case[11] ]],
                          const.var    = options_const.var[[    case[12] ]],
                          inits        = options_inits[[        case[13] ]],
                          n.iter       = options_n.iter[[       case[14] ]],
                          n.chains     = options_n.chains[[     case[15] ]],
                          n.burn       = options_n.burn[[       case[16] ]],
                          cachesize    = options_cachesize[[    case[17] ]],
                          quiet        = options_quiet[[        case[18] ]]),
                  throws_error(msg))
  }
}

#Table containing all of the test cases
#if(FALSE){
#  testcases <- read.csv("tests/testthat/testcases.csv")
  testcases <- read.csv("testcases.csv")
  for(index in 2:nrow(testcases)){
    case <- as.numeric(testcases[index,-c(1,ncol(testcases))])
    test_that(
      paste0("Scenario: ", testcases[index,1], " --- ", testcases[index,ncol(testcases)]),
      PeriodCPT_TEST(case)
    )
  }
#}

##TODO
#check iter input
#  -- periodlength
#  -- minseglen
#  -- n.chains
#  -- DOTDOTDOT!!!
#  -- (internal sim) Mprior
#  -- (internal sim) Mhyp
#  -- (internal sim) spread
