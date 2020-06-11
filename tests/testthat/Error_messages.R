##List of Error Messages!!!
###########################

ErrorMessages <- c(
  "'arg' should be one of “bern”, “pois”, “norm”, “mean”, “var”",
  "'arg' should be one of “pois”, “unif”",
  "Data is missing.",
  "Data must not contain NA missing values.",
  "Period length is not defined either via data as `ts` object or explicitly given as input.",
  "Mhyp specified incorrectly for Mprior `$(sub_st)$options_Mprior[[case[5]]]$(sub_ed)$`.",
  "Implementation Error: Mprior `$(sub_st)$options_Mprior[[case[5]]]$(sub_ed)$` is not supported.",     ##<= is this a printed error!!!
  "Hyper-parameter `spread` specified incorrectly.",
  "MCMC option - n.iter not specified.",
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


###extend
#"Unexpected class of `object`."
#"Invalid newiters."
#"newiters must be a positive integer."
#"Cannot determine last sample to initiate next batch."

##summarise
"Chains have been summarised, cannout use summarised data to determine last chain samples."
