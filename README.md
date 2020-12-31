# PLSR modeling for the estimation of plant functional traits
This repository contains example scripts illustrating "best-practices" for fitting, evaluating, and reporting leaf-level spectra-trait PLSR models. These scripts encompass several possibilities that you may encounter doing PLSR. Start by reading *Burnett et al. 2020*, then work through the scripts or vignettes.

### Article citation:
Burnett AC, Anderson J, Davidson KD, Ely KS, Lamour J, Li Q, Morrison BD, Yang D, Rogers A, Serbin SP (2020) A best-practice guide to predicting plant traits from leaf-level hyperspectral data using partial least squares regression. Journal of Experimental Botany. In Review.

### Source code citation:
[![DOI](https://zenodo.org/badge/222699149.svg)](https://zenodo.org/badge/latestdoi/222699149)

### EcoSML
https://ecosml.org/package/github/TESTgroup-BNL/PLSR_for_plant_trait_prediction

### Script authors:
Julien Lamour, Jeremiah Anderson, Ken Davidson, Shawn P. Serbin

### Depends: 
ggplot2 (>= 3.3.2), remotes (>= 2.2.0), devtools (>= 2.3.1), readr (>= 1.3.1), RCurl (>= 1.98-1.2), httr (>= 1.4.2), pls (>= 2.7-2), dplyr (>= 1.0.1), reshape2 (>= 1.4.4), here (>= 0.1), plotrix (>= 3.7-8), gridExtra (>= 2.3), scales(>= 1.1.1)

## Contains:
1. Example R script files that illustrate the "best-practices" of PLSR model fitting for the estimation of leaf functional traits with reflectance spectroscopy
    * _spectra-trait_kit_sla_plsr_example.R_ Small dataset looking at SLA with some data cleaning (removal of NA's and suspect high values)
    * _spectra-trait_neon_lma_plsr_example.R_ Large dataset looking at LMA with multiple grouping variables, very slow (>6000 observations)
    * _spectra-trait_reseco_leafN_plsr_example.R_ Small dataset looking at leaf nitrogen content
    * _spectra-trait_reseco_lma_plsr_example.R_ Small dataset looking at LMA
    * _simple_spectra-trait_plsr_example.R_ Basic PLSR example using a large dataset
    * _pull_data_from_ecosis_ Quick example of how to pull data from EcoSIS and plot it

2. Non-CRAN or external library R functions used in the example PLSR model fitting scripts provided in the "functions.R" file
    * _get_ecosis_data()_ Function to pull data from the EcoSIS database (ecosis.org) using their application programmer interface (API)
    * _create_data_split()_ Randomly splits data into calibration and validation datasets based on grouping variables.  'base' option is slow but verbose.  'dplyr' is fast and quiet.
    * _f.plot.spec()_ Function to generate spectral plot with mean, min/max and 95% confidence intervals
    * _find_optimal_components()_ Finds optimum number of components for PLSR.  'pls' chooses the model with fewest components that is still less than one standard error away from the overall best model. 'first plateau' chooses the first component that gives statistically (t-test) the same result as the following component.  'firstMin' finds the first component that gives statistically (t-test) the same result as the overall best model.
    * _pls_permutation()_ Generates PLSR model permutation analysis ensembles for opimal component selection and uncertainty analysis.  Currently called by _find_optimal_components()_
    * _f.plot.coef()_ Plots PLSR model coefficients with uncertainty envelope
    * _f.coef.valid()_ Returns the intercept and the coefficients of the jackknife permutation analysis.
  
3. Example Rmarkdown vignettes illustrating the various PLSR model fitting examples

### Linked dataset citations, DOIs, and EcoSIS IDs/URLs: <br>
1) Leaf reflectance plant functional gradient IFGG/KIT <br>
Target variable: SLA <br>
EcoSIS URL: https://ecosis.org/package/leaf-reflectance-plant-functional-gradient-ifgg-kit <br>
EcoSIS ID: 3cf6b27e-d80e-4bc7-b214-c95506e46daa <br>
Rpubs example output: https://rpubs.com/sserbin/705140

2) Fresh Leaf Spectra to Estimate LMA over NEON domains in eastern United States <br>
Target variable: LMA <br>
EcoSIS URL: https://ecosis.org/package/fresh-leaf-spectra-to-estimate-lma-over-neon-domains-in-eastern-united-states <br>
EcoSIS ID: 5617da17-c925-49fb-b395-45a51291bd2d <br>
DOI: https://doi.org/doi:10.21232/9831-rq60 <br>
Rpubs example output: https://rpubs.com/sserbin/665529

3) Leaf spectra of 36 species growing in Rosa rugosa invaded coastal grassland communities in Belgium <br>
Target variable: LMA, leaf nitrogen <br>
EcoSIS URL: https://ecosis.org/package/leaf-spectra-of-36-species-growing-in-rosa-rugosa-invaded-coastal-grassland-communities-in-belgium <br>
EcoSIS ID: 9db4c5a2-7eac-4e1e-8859-009233648e89 <br>
DOI: https://doi.org/doi:10.21232/9nr6-sq54 <br>
Rpubs LMA example output: https://rpubs.com/sserbin/665512 <br>
Rpubs LeafN example output: https://rpubs.com/sserbin/665516


