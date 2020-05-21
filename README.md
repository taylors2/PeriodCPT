## Package information

This package implements the within period changepoint detection algorith that is based on a paper that is under review which focuses on periodic (eg daily) binary time series data.

The PeriodCPT package is under development with the aim of including 
utility for different data types and appropriate changes in different structures. This package is not yet published on CRAN and will not be until a satisfactory version has been made. In the meantime, use with great causion and please report any bugs.

## List of task for PeriodCPT
- [x] Put package on GitHub.
- [x] Remove utilities not to be included in final package.
- [x] Audit functions for which cases are suppoted.
- [ ] Add the following missing/extra functionality.
  - [x] Poisson, Fit Measure
  - [x] Bernoulli, Fit Measure
  - [x] Normal mean, Fit Measure
  - [x] Normal var, Fit Measure
  - [ ] Generic functions.
- [x] Correct distribution name from Binomial to Bernoulli
  - [ ] Add functionality for Binomial (?) [some code is there, full implementation not yet avilable]
- [x] Deal with various input data classes (eg. ts objects)
  - [x] CircCT() - edit time input to specify period length (Default NULL, take time from ts object).
- [x] Make C code generic for the addition of new distribution options.
- [x] Update manual pages.
- [ ] Create tests.
- [ ] Plot function - adapt period/time axis.
- [x] Warning/Error messages.
- [ ] Use position in period index rather than based on 24hr.
- [ ] Convergence checks. (Related to MSL) 
- [ ] Speed test.
- [ ] Allow data being specifed as zoo object (nb, may need coredata() to extract time series data.)
- [ ] plot functions
- [ ] quantile functions
- [ ] Bug Fixes!!!
NOTE: Do not edit 'MODE' item in CircCPT output.
