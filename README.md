## CHECKS

[![Build Status](https://travis-ci.com/taylors2/PeriodCPT.svg?token=m8aKb1xsE6ZMGUKbKdHG&branch=master)](https://travis-ci.com/taylors2/PeriodCPT) [![Codecov test coverage](https://codecov.io/gh/taylors2/PeriodCPT/branch/master/graph/badge.svg)](https://codecov.io/gh/taylors2/PeriodCPT)

## Package information

This package implements the within period changepoint detection algorith that is based on a paper that is under review which focuses on periodic (eg daily) binary time series data.

The PeriodCPT package is under development with the aim of including 
utility for different data types and appropriate changes in different structures. This package is not yet published on CRAN and will not be until a satisfactory version has been made. In the meantime, use with great causion and please report any bugs.

## List of task for PeriodCPT
- [ ] Create tests.
  - [ ] PeriodCPT()
  - [ ] PeriodCPT.extend()
  - [ ] summarise
  - [ ] plot
  - [ ] quantiles
- [ ] plot functions
- [ ] quantile functions
- [ ] Convergence checks. (Related to MSL) 
- [ ] Add functionality for Binomial (?) [some code is there, full implementation not yet avilable]
- [ ] Allow data being specifed as zoo object (nb, may need coredata() to extract time series data.)
