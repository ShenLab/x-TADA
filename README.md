# x-TADA
x-TADA is a new model to integrate functional genomics data in identifying de novo risk genes. Inspired from TADA (He et al. 2013) and extTADA (Nguyen et al. 2017) but different. In additional to the observation of damage variants, it takes the expression level of genes during embryo development as another orthogonal obervation. 

  We hypothesized that disease risk genes should:
  
    1) Harbor enriched damage de novo variants in probands compared to non-affected population. 
    2) Show high expression during corresponding organ development.
    
  The model look like this:
  
  <img src="https://github.com/ShenLab/x-TADA/blob/master/x-TADA.model.png?raw=true" width="512">
  
ϕ_0 and ϕ_1 are hyperparameters to estimate. N is the observed damage variant number in patient cohorts. S is the expression level for each gene, in this CHD practice it is mouse E14.5 heart expression rank percentile.

## Environment:
  
  `R>=3.6.0`
  
  `packages: stan, dplyr, funr`

## Usage:

  Source all the R files in `model/` directory.
  
  Modify the `demo.R` with your input. See input format in demo.R
  
  Run `source(demo.R)`.
  
## Result:
  
  A list object with three attributes:
  
    1) mcmcDD
  
      mcmcDD raw result
  
    2) pars0
  
      statistical summary of mcmcDD result
  
    3) dataFDR
  
      posterior probability, qvalues, etc for each gene.
  
