Environment:
  
  R>=3.6.0
  packages: stan, dplyr 

Usage:

  Source all the R files in model/ directory.
  
  Modify the demo.R with your input. See input format in demo.R
  
  Source demo.R.
  
Result:
  
  A list object with three attributes:
  
  1) mcmcDD
  
  mcmcDD raw result
  
  2) pars0
  
  statistical summary of mcmcDD result
  
  3) dataFDR
  
  posterior probability, qvalues, etc for each gene.
  
