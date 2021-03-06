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
    data <- ts(data = rnorm(cycles*periodlength, mean = 0, sd = 2),
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
    data <- ts(data = rnorm(cycles*periodlength, mean = 5, sd = 0.5),
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
  }else if(caseid == 15){
    data <- ts(matrix(data,nc=1), frequency = frequency(data))
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
                           NA, NULL, -1, 0.1, 1+24, "A", 1) #invalid
options_Mprior <- list("pois","unif",          #valid
                       "invalid", 1, NA, NULL) #invalid
options_Mhyp <- list(1, #valid
                     NA, NULL, 0, -1, "A") #invalid    (nb: only applicable when Mprior == "pois")
options_spread <- list(1, #valid
                       NA, NULL, 0, -1, "A") #invalid
options_param.a <- options_param.b <- options_param.m <-
  options_param.c <- options_const.var <- list(1, #valid
                                           NA, NULL, 0, -1, "A") #invalid
inits_fn1   <- function( ... ){
  return(1)
}

inits_fn2   <- function(pcpt.object, chain.index, ... ){
  if(chain.index == 1){
      return(c(1, round(periodlength(object = pcpt.object)/2)))
  }else{
    return(1)
  }
}

inits_fn3 <- function(unspecified_input){
  return(unspecified_input)
}

inits_fn4 <- function(unspecified_input, ...){
  return(unspecified_input)
}

options_inits <- list(
  NULL,                            ##should simulate from prior
  "ends",                          ##should over-ride n.chain
  list(1,1),                       #Valid if n.chain <= 2
  c(1,2),                          #valid if minseglen == 1
  c(0,24/2),                       #invalid
  c(1,24+1),                       #invalid
  c(24,24/2),                      #invalid
  rep(1,24+1),                     #invalid
  c(1, 24/2 + 0.5),                #invalid
  c("1", "12"),                    #invalid
  list(1,numeric(0)),              #invalid
  c(1,NA,24/2),                    #invalid
  inits_fn1,                 ##should define value from function
  inits_fn2,                 ##should define value from function
  inits_fn3,                  #invalid
  list(1),                     #invalid for nchains>1
  list(1,"1"),
  inits_fn4                  #invalid
)

#options_DOTDOTDOT <- list(inits_fn_args = 0,  #valid
#                          pcpt.object=new("pcpt"), chain.index=1) #invalid

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
  "gettextf(\"\'arg\' should be one of %s\", paste(dQuote(c('bern', 'pois', 'norm', 'mean', 'var')), collapse = ', '))",
  "gettextf(\"\'arg\' should be one of %s\", paste(dQuote(c('pois','unif')), collapse = ', '))",
  "Data is missing.",
  "Data must not contain NA missing values.",
  "Period length is not defined either via data as `ts` object or explicitly given as input.",
  "Mhyp specified incorrectly for Mprior `$(sub_st)$options_Mprior[[case[5]]]$(sub_ed)$`.",
  "Unexpected class of `object`.",
  "Hyper-parameter `spread` specified incorrectly.",
  "MCMC option - n.iter specified incorrectly.",
  "MCMC option - n.chains specified incorrectly.",    ###10
  "MCMC option - n.burn specified incorrectly.",
  "MCMC option - cachesize specified incorrectly.",
  "MCMC option - quiet specified incorrectly.",
  "Argument `newiters` must be a single positive integer.",
  "Cannot pass arguments 'pcpt.object' and 'chain.index' from '...' into inits(). These inputs to inits() are managed internally.",
  "Argument `object` does not appear to be an output from a PeriodCPT function.",
  "Incorrect number of initial values for specified number of chains.",
  "Class of at least one inits is not numeric.",
  "Incorrect number of within period changepoints specified by inits.",
  "In inits, within period cpts must be whole numbers.",    ###20
  "In inits, within period cpts must be between 1 and period length.",
  "In inits, within period cpts are not ordered or do not satisfy minimum segment length condition.",
  "Hyper-parameter `param.a` specified incorrectly.",
  "Hyper-parameter `param.b` specified incorrectly.",
  "Hyper-parameter `param.c` specified incorrectly.",
  "Hyper-parameter `param.m` specified incorrectly.",
  "Known constant `const.var` is specified incorrectly.",
  "Data is invalid for '$(sub_st)$options_distribution[[case[1]]]$(sub_ed)$' sampling distribution.",
  "Argument `all` is not a single logical value.",
  "argument \"object\" is missing, with no default",         ###30
  "Cannot access information for index `$(sub_st)$index$(sub_ed)$`. It is likely that the contents accoss multiple chains has been combined into the first index.",
  "Period length specified incorrectly.",
  "Minimum segment length specifed incorrectly.",
  "Argument 'index' missing with no default.",
  "Argument 'index' must be a single numeric value.",
  "Index `$(sub_st)$options_index[[i]]$(sub_ed)$` not found in results list.",
  "Probabilities must be a numerical vector between 0 and 1.",
  "",
  "Can only assign logical to summarised slot.",
  "Assignment to nsegparam slot is not a single positive integer.",    ###40
  "Length of data is too short or period length is too long.",
  "Period length must be greater than 1.",
  "Minimum segment length longer than period length.",
  "'arg' must be NULL or a character vector",
  "Mprior cannot be NULL.",
  "MCMC option - n.iter not specified.",
  "Unexpected class of `object`.",
  "Inits must not contain missing NA values.",
  "Inital values have been specified incorrectly.",
  "Passing unused arguments to inits function, consider adding '...' to the inputs of your function.",  ###50
  "\n  argument \"unspecified_input\" is missing, with no default\n"
  )

##Function to perform the testthat commands
PeriodCPT_TEST <- function(case){
  ##Generate data
  data <- make_test_data(options_distribution[[case[1]]], case[2])

  FNSTR <- paste0("PeriodCPT(data = data, ",
            "distribution = options_distribution[[ case[01] ]], ",
            "periodlength = options_periodlength[[ case[03] ]], ",
            "minseglen    = options_minseglen[[    case[04] ]], ",
            "Mprior       = options_Mprior[[       case[05] ]], ",
            "Mhyp         = options_Mhyp[[         case[06] ]], ",
            "spread       = options_spread[[       case[07] ]], ",
            "param.a      = options_param.a[[      case[08] ]], ",
            "param.b      = options_param.b[[      case[09] ]], ",
            "param.c      = options_param.c[[      case[10] ]], ",
            "param.m      = options_param.m[[      case[11] ]], ",
            "const.var    = options_const.var[[    case[12] ]], ",
            "inits        = options_inits[[        case[13] ]], ",
            "n.iter       = options_n.iter[[       case[14] ]], ",
            "n.chains     = options_n.chains[[     case[15] ]], ",
            "n.burn       = options_n.burn[[       case[16] ]], ",
            "cachesize    = options_cachesize[[    case[17] ]], ")
  if(case[length(case)-1] == 2){
    FNSTR <- paste0(FNSTR," pcpt.object = new(\"pcpt\"), ")
  }else if(case[length(case)-1] == 3){
    FNSTR <- paste0(FNSTR, "  chain.index = 1, ")
  }else if(case[length(case)-1] == 4){
    FNSTR <- paste0(FNSTR, "  unspecified_input = 1, ")
  }
  FNSTR <- paste0(FNSTR, "quiet        = options_quiet[[        case[18] ]])")


  if(case[length(case)]==0){
      #Expect function to be valid, returning a "pcpt" class object
      expect_s4_class(eval(parse(text = FNSTR)), "pcpt")
  }else{
      #Expect test to return error message
      msg <- ErrorMessages[case[length(case)]]
      if(case[length(case)] <= 2){
        msg <- eval(parse(text = msg))
      }else if( case[length(case)] != 51){
        msg <- gsub('$(sub_st)$', '",', msg, fixed=TRUE)
        msg <- gsub('$(sub_ed)$', ',"', msg, fixed=TRUE)
        msg <- eval(parse(text = paste0('paste0("',msg,'")')))
      }
      expect_that(eval(parse(text = FNSTR)), throws_error(msg, fixed = TRUE))
  }
}

options_summarise_slot <- list(TRUE, FALSE, c(TRUE, FALSE), NA, NULL, 0, -1, "A", 0.5)
summarise_test_function <- function(value){
  x <- new("pcpt")
  summarised(x) <- value
  return(x)
}

options_nsegparam_slot <- list(1, TRUE, NA, NULL, 0, -1, "A", 1.5, c(1, 2))
nsegparam_test_function <- function(value){
  x <- new("pcpt")
  PeriodCPT:::nsegparam(x) <- value
  return(x)
}

result_index_function <- function(LIST, index){
  x <- new("pcpt")
  results(x) <- LIST
  result(x, index) <- 1
  return(x)
}
RES_LIST <- RES_LIST_ERR <- list(0,0)
names(RES_LIST) <- as.character(1:2)

options_index <- list(0,-1,0.5,"A",NA,NULL,c(1,2),c(2,3))


################################################################
################################################################

####TESTTHAT: PeriodCPT()
RUN <- TRUE
#RUN <- FALSE
if(RUN){
  testcases <- read.csv("testcases.csv")
}else{
  library(testthat)
  testcases <- read.csv("tests/testthat/testcases.csv")
}

test_that("Minimal example", {expect_s4_class(PeriodCPT(data = make_test_data("bern", 1), distribution = "bern"), "pcpt")})
test_that("Bad class", {expect_that(PeriodCPT:::PeriodCPT.main(1), throws_error(ErrorMessages[47]))})
test_that("Bad data - master", {expect_that(PeriodCPT(),      throws_error(ErrorMessages[3]))})
test_that("Bad data - bern",   {expect_that(PeriodCPT.bern(), throws_error(ErrorMessages[3]))})
test_that("Bad data - mean",   {expect_that(PeriodCPT.mean(), throws_error(ErrorMessages[3]))})
test_that("Bad data - norm",   {expect_that(PeriodCPT.norm(), throws_error(ErrorMessages[3]))})
test_that("Bad data - pois",   {expect_that(PeriodCPT.pois(), throws_error(ErrorMessages[3]))})
test_that("Bad data - var",    {expect_that(PeriodCPT.var(),  throws_error(ErrorMessages[3]))})

for(index in 1:nrow(testcases)){
  if(testcases[index,1] > 0){
    case <- as.numeric(testcases[index,-c(1,ncol(testcases))])
    test_that(
      paste0("Scenario: ", testcases[index,1], " --- ", testcases[index,ncol(testcases)]),
      PeriodCPT_TEST(case)
    )
  }
}

for(i in 1:length(options_summarise_slot)){
  test_that(paste0("summarise_slot: ",i), {
    if(i<=3){
      expect_s4_class(summarise_test_function(options_summarise_slot[[i]]), "pcpt")
    }else{
      expect_that(summarise_test_function(options_summarise_slot[[i]]),
                  throws_error(ErrorMessages[39], fixed = TRUE))
    }
  })
}

for(i in 1:length(options_nsegparam_slot)){
  test_that(paste0("nsegparam_slot: ",i), {
    if(i==1){
      expect_s4_class(nsegparam_test_function(options_nsegparam_slot[[i]]), "pcpt")
    }else{
      expect_that(nsegparam_test_function(options_nsegparam_slot[[i]]),
                  throws_error(ErrorMessages[40], fixed = TRUE))
    }
  })
}

test_that("Results_list: 1",expect_s4_class(result_index_function(RES_LIST,1),"pcpt"))
test_that("Results_list: 2",expect_s4_class(result_index_function(RES_LIST,"1"),"pcpt"))
test_that("Results_list: 3",expect_that(result_index_function(RES_LIST,"A"),    throws_error("Index `A` not found in results list.")))
test_that("Results_list: 4",expect_that(result_index_function(RES_LIST_ERR,1),  throws_error("Index `1` not found in results list.")))
test_that("Results_list: 5",expect_that(result_index_function(RES_LIST_ERR,"1"),throws_error("Index `1` not found in results list.")))
test_that("Results_list: 6",expect_that(result_index_function(RES_LIST_ERR,"A"),throws_error("Index `A` not found in results list.")))

test_that("Edits with periodlength and minseglen", expect_that({
    x <- new("pcpt")
    data.set(x) <- ts(rbinom(100, size = 1, prob = 0.5), frequency = 5)
    periodlength(x) <- 20
    minseglen(x) <- 7 #defines npcpt.max
    periodlength(x) <- NULL
  }, throws_error( ErrorMessages[43] )))


test_that("Edits with periodlength and minseglen v2", expect_equal({
  tmp <- c(NA,NA)
  x <- new("pcpt")
  data.set(x) <- ts(rbinom(100, size = 1, prob = 0.5), frequency = 20)
  minseglen(x) <- 2 #defines npcpt.max
  tmp[1] <- npcpts.max(x)
  periodlength(x) <- 10 #adjusts npcpt.max
  tmp[2] <- npcpts.max(x)
  tmp
}, c(10, 5)))



#######################
x <- ts(sample(c(0,1), size = 240, replace = TRUE), frequency = 24)
ans <- PeriodCPT(data = x, distribution = "bern", n.iter = 100, quiet = TRUE, inits = "ends")
test_that("Accessing Results - missing index",
          expect_that(result(ans),throws_error(ErrorMessages[34])))
for(index in 1:n.chains(ans)){
  test_that(paste0("Accessing Results - good, chain ",index),
            expect_is(result(ans, index),"matrix"))
}
for(i in 1:length(options_index)){
  if(i<=3){
    msg <- ErrorMessages[36]
    msg <- gsub('$(sub_st)$', '",', msg, fixed=TRUE)
    msg <- gsub('$(sub_ed)$', ',"', msg, fixed=TRUE)
    msg <- eval(parse(text = paste0('paste0("',msg,'")')))
  }else{
    msg <- ErrorMessages[35]
  }
  test_that(paste0("Accessing Results - bad index (", i,")"),
            expect_that(result(ans, options_index[[i]]),throws_error(msg,fixed = TRUE)))
}

test_that("Accessing results", expect_is(results(ans),"list"))


####TESTTHAT: PeriodCPT_extend()

x <- ts(sample(c(0,1), size = 240, replace = TRUE), frequency = 24)
ans <- PeriodCPT(data = x, distribution = "bern", n.iter = 100, quiet = TRUE)

test_that("Extend - minimal use",
          expect_s4_class({PeriodCPT_extend(object = ans)},"pcpt"))
test_that("Extend - specify valid newiters",
          expect_s4_class({PeriodCPT_extend(object = ans, newiters = options_n.iter[[1]])},
                         "pcpt"))
test_that("Extend - valid with summarise",
          expect_s4_class({PeriodCPT_extend(object = summarise_chains(ans),
                                           newiters = options_n.iter[[1]])},"pcpt"))

test_that("Extend - missing input",
         expect_that({PeriodCPT_extend(newiters = options_n.iter[[1]])},
                     throws_error(ErrorMessages[30])))
test_that("Extend - class error",
         expect_that({PeriodCPT_extend(object = 1, newiters = options_n.iter[[1]])},
                     throws_error(ErrorMessages[7])))
for(i in 2:length(options_n.iter)){
  test_that(paste0("Extend - invalid newiters (",i,")"),
    expect_that({PeriodCPT_extend(ans, newiters = options_n.iter[[i]])},
                throws_error(ErrorMessages[14])))
}
test_that("Extend - invalid pcpt object",
  expect_that({PeriodCPT_extend(object = new("pcpt"))}, throws_error(ErrorMessages[16])))

####TESTTHAT: summarise_chains()

options_all <- list(TRUE, FALSE,  #valid
                     NA,NULL,0,-1,"A",0.5,c(TRUE,FALSE)) #invalid
x <- ts(sample(c(0,1), size = 240, replace = TRUE), frequency = 24)
ans <- PeriodCPT(data = x, distribution = "bern", n.iter = 100, quiet = TRUE)

test_that("SummeriSe - missing input",
          expect_that(summarise_chains(all=TRUE), throws_error(ErrorMessages[30])))
test_that("SummeriSe - class error",
          expect_that(summarise_chains(object = 1), throws_error(ErrorMessages[7])))
test_that("SummeriSe - object not an output from PeriodCPT",
          expect_that(summarise_chains(object = new("pcpt")), throws_error(ErrorMessages[16])))
test_that("SummeriZe - missing input",
          expect_that(summarize_chains(), throws_error(ErrorMessages[30])))
test_that("SummeriZe - class error",
          expect_that(summarize_chains(object = 1), throws_error(ErrorMessages[7])))
test_that("SummeriZe - object not an output from PeriodCPT",
          expect_that(summarize_chains(object = new("pcpt")), throws_error(ErrorMessages[16])))

for(i in 1:length(options_all)){
  if(i<=2){
  test_that("SummariSe -- one chain", expect_s4_class(
    summarise_chains(ans, all = options_all[[i]]), "pcpt"))
  test_that("SummariZe -- one chain", expect_s4_class(
      summarize_chains(ans, all = options_all[[i]]), "pcpt"))
  }else{
    test_that(paste0("SummariSe -- one chain, invalid `all` argument (",i,")"),
      expect_that(summarise_chains(ans, all = options_all[[i]]),
                  throws_error(ErrorMessages[29])))
    test_that(paste0("SummariZe -- one chain, invalid `all` argument (",i,")"),
      expect_that(summarize_chains(ans, all = options_all[[i]]),
                  throws_error(ErrorMessages[29])))
  }
}

x <- ts(sample(c(0,1), size = 240, replace = TRUE), frequency = 24)
ans <- PeriodCPT(data = x, distribution = "bern", n.iter = 100, n.chains = 2, quiet = TRUE)

test_that("SummariSe -- two chains", expect_s4_class(summarise_chains(ans), "pcpt"))
test_that("SummariZe -- two chains", expect_s4_class(summarize_chains(ans), "pcpt"))
test_that("SummariSe -- two chains, twice",
          expect_s4_class({
            tmp <- summarise_chains(ans)
            summarise_chains(tmp)
          },"pcpt"))
test_that("SummariSe -- two chains, individually and then combine",
          expect_s4_class({
            tmp <- summarise_chains(ans, all = FALSE)
            summarise_chains(tmp, all = TRUE)
          },"pcpt"))
test_that("SummariZe -- two chains, twice",
          expect_s4_class({
            tmp <- summarize_chains(ans)
            summarize_chains(ans)
          },"pcpt"))
test_that("SummariZe -- two chains, individually and then combine",
          expect_s4_class({
            tmp <- summarize_chains(ans, all = FALSE)
            summarize_chains(tmp, all = TRUE)
          },"pcpt"))

x <- ts(sample(c(0,1), size = 240, replace = TRUE), frequency = 24)
ans <- PeriodCPT(data = x, distribution = "bern", n.iter = 100, quiet = TRUE, inits = "ends")
test_that("SummariSe -- two chains, ncol(tab1) < ncol(tab2)",
          expect_s4_class(summarise_chains(ans, all = TRUE),"pcpt"))
test_that("SummariZe -- two chains, ncol(tab1) < ncol(tab2)",
          expect_s4_class(summarize_chains(ans, all = TRUE),"pcpt"))
inits_rev <- list(1:periodlength(ans),1)
ans <- PeriodCPT(data = x, distribution = "bern", n.iter = 100, quiet = TRUE, inits = inits_rev, n.chains = 2)
test_that("SummariSe -- two chains, ncol(tab1) > ncol(tab2)",
          expect_s4_class(summarise_chains(ans, all = TRUE),"pcpt"))
test_that("SummariZe -- two chains, ncol(tab1) > ncol(tab2)",
          expect_s4_class(summarize_chains(ans, all = TRUE),"pcpt"))

x <- ts(sample(c(0,1), size = 240, replace = TRUE), frequency = 24)
ans <- PeriodCPT(data = x, distribution = "bern", n.iter = 100, quiet = TRUE, inits = "ends")
ans <- summarise_chains(ans, all = TRUE)
for(index in 1:n.chains(ans)){
  if(index==1){
    test_that(paste0("SummariS/Ze -- two chains, accessing results (",index,")."),
          expect_is(result(ans,index),"matrix"))
  }else{
    msg <- ErrorMessages[31]
    msg <- gsub('$(sub_st)$', '",', msg, fixed=TRUE)
    msg <- gsub('$(sub_ed)$', ',"', msg, fixed=TRUE)
    msg <- eval(parse(text = paste0('paste0("',msg,'")')))
    test_that(paste0("SummariS/Ze -- two chains, accessing results (",index,")."),
              expect_that(result(ans,index),throws_error(msg, fixed = TRUE)))
  }
}

for(dist in options_distribution){
  if(dist == "invalid") next
  data <- make_test_data(dist, 1)
  ans <- PeriodCPT(data, dist, quiet = TRUE, n.iter = 100)
  test_that(paste0("SummariSe - ",dist), expect_s4_class(summarise_chains(ans),"pcpt"))
  test_that(paste0("SummariZe - ",dist), expect_s4_class(summarise_chains(ans),"pcpt"))
}


####TESTTHAT: table_npcpt()

x <- ts(sample(c(0,1), size = 240, replace = TRUE), frequency = 24)
ans <- PeriodCPT(data = x, distribution = "bern", n.iter = 100, quiet = TRUE, inits = "ends")
ans2 <- summarise_chains(ans, all = FALSE)
ans3 <- summarise_chains(ans, all = TRUE)
test_that("TableM - missing input",
          expect_that(table_npcpt(), throws_error(ErrorMessages[30])))
test_that("TableM - class error",
          expect_that(table_npcpt(1), throws_error(ErrorMessages[7])))
test_that("TableM - object not an output from PeriodCPT",
          expect_that(table_npcpt(object = new("pcpt")), throws_error(ErrorMessages[16])))
test_that("TableM - Good, using raw chains",
          expect_is(table_npcpt(object = ans), "table"))
test_that("TableM - Good, using separtate summarised chains",
          expect_is(table_npcpt(object = ans2), "table"))
test_that("TableM - Good, using overall summarised chains",
          expect_is(table_npcpt(object = ans3), "table"))


####TESTTHAT: table_pcpt()

x <- ts(sample(c(0,1), size = 240, replace = TRUE), frequency = 24)
ans <- PeriodCPT(data = x, distribution = "bern", n.iter = 100, quiet = TRUE, inits = "ends")
ans2 <- summarise_chains(ans, all = FALSE)
ans3 <- summarise_chains(ans, all = TRUE)
test_that("TableTAU - missing input",
          expect_that(table_pcpt(), throws_error(ErrorMessages[30])))
test_that("TableTAU - class error",
          expect_that(table_pcpt(1), throws_error(ErrorMessages[7])))
test_that("TableTAU - object not an output from PeriodCPT",
          expect_that(table_pcpt(object = new("pcpt")), throws_error(ErrorMessages[16])))
test_that("TableTAU - Good, using raw chains",
          expect_is(table_pcpt(object = ans), "table"))
test_that("TableTAU - Good, using separtate summarised chains",
          expect_is(table_pcpt(object = ans2), "table"))
test_that("TableTAU - Good, using overall summarised chains",
          expect_is(table_pcpt(object = ans3), "table"))





test_that(paste0("Quantile - ",dist,", not PeriodCPT output"),
          expect_that(quantile(new("pcpt")),throws_error(ErrorMessages[16])))
for(dist in unlist(options_distribution)){
  if(dist == "invalid") next
  ans <- PeriodCPT(make_test_data(dist, 1), distribution = dist, quiet = TRUE, n.iter = 1000)
  test_that(paste0("Quantile - ",dist,", default, not summarised"),
            expect_is(quantile(ans),"matrix"))
  ans <- summarise_chains(ans, all = TRUE)
  test_that(paste0("Quantile - ",dist,", default, summarised"),
            expect_is(quantile(ans),"matrix"))
  test_that(paste0("Quantile - ",dist,", default, other ignored imput"),
            expect_is(quantile(ans, junk = "RUBBISH"),"matrix"))
  test_that(paste0("Quantile - ",dist,", p = 0.5, summarised"),
            expect_is(quantile(ans, probs = 0.5),"matrix"))
  test_that(paste0("Quantile - ",dist,", p = 95%CI, summarised"),
            expect_is(quantile(ans, probs = c(0.025,0.975)),"matrix"))

  test_that(paste0("Quantile - ",dist,", p = NA, summarised"),
            expect_that(quantile(ans, probs = NA),throws_error(ErrorMessages[37])))
  test_that(paste0("Quantile - ",dist,", p = NULL, summarised"),
            expect_that(quantile(ans, probs = NULL),throws_error(ErrorMessages[37])))
  test_that(paste0("Quantile - ",dist,", p = too big, summarised"),
            expect_that(quantile(ans, probs = c(0.025,1.1)),throws_error(ErrorMessages[37])))
  test_that(paste0("Quantile - ",dist,", p = too small, summarised"),
            expect_that(quantile(ans, probs = c(-0.1, 0.975)),throws_error(ErrorMessages[37])))
  test_that(paste0("Quantile - ",dist,", p = char, summarised"),
            expect_that(quantile(ans, probs = c("0.025", "0.975")),throws_error(ErrorMessages[37])))

  ans <- PeriodCPT(make_test_data(dist, 2), distribution = dist, quiet = TRUE, n.iter = 1000)
  test_that(paste0("Quantile - ",dist,", default, not summarised, single segment"),
            expect_is(quantile(ans),"matrix"))
}

