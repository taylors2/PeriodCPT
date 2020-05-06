## Package information

This package implements the within period changepoint detection algorith that is based on a paper that is under review which focuses on periodic (eg daily) binary time series data.

The PeriodCPT package is under development with the aim of including 
utility for different data types and appropriate changes in different structures. This package is not yet published on CRAN and will not be until a satisfactory version has been made. In the meantime, use with great causion and please report any bugs.

## List of task for PeriodCPT
- [x] Put package on GitHub.
- [x] Remove utilities not to be included in final package.
- [ ] Audit functions for which cases are suppoted.
- [ ] Add the following missing/extra functionality.
  - [ ] Poisson, Fit Measure
  - [ ] Bernoulli, Fit Measure
  - [ ] Normal mean, Fit Measure
  - [ ] Normal var, Fit Measure
  - [ ] Generic functions.
- [ ] Correct distribution name from Binomial to Bernoulli
  - [ ] Add functionality for Binomial (?)
- [ ] Deal with various input data classes (eg. ts objects)
  - [ ] CircCT() - edit time input to specify period length (Default NULL, take time from ts object).
- [ ] Make C code generic for the addition of new distribution options.
- [ ] Update manual pages.
- [ ] Create tests.
- [ ] Plot function - adapt period/time axis.
- [ ] Warning/Error messages.
- [ ] Use position in period index rather than based on 24hr.
- [ ] Convergence checks. (Related to MSL)
- [ ] Speed test.
- [ ] Allow data being specifed as zoo object (nb, may need coredata() to extract time series data.)
NOTE: Do not edit 'MODE' item in CircCPT output.
