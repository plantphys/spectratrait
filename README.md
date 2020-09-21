# PLSR modeling for the estimation of plant functional traits
Tutorial and example scripts illustrating how to fit, evaluate, and report spectra-trait reflectance spectroscopy PLSR models


### Contains:
1) Example R script files that illustrate the "best practices" of PLSR model fitting for the estimation of leaf functional traits with reflectance spectroscopy
  + _expanded_spectra-trait_kit_lma_plsr_example.R_  Small dataset looking at LMA with some data cleaning (removal of NA's and suspect high values)
  + _expanded_spectra-trait_neon_lma_plsr_example.R_ Large dataset looking at LMA with multiple grouping variables, very slow (>6000 observations)
  + _expanded_spectra-trait_reseco_leafN_plsr_example.R_ Small dataset looking at leaf nitrogen
  + _expanded_spectra-trait_reseco_lma_plsr_example.R_ Small dataset looking at LMA
  + _simple_spectra-trait_plsr_example.R_ Basic PLSR example using a large dataset
  + _pull_data_from_ecosis_ Quick example of how to pull data from EcoSIS and plot it

2) R functions used in the example PLSR model fitting scripts
  + _create_data_split()_ Randomly plits data into calibration and validation datasets based on grouping variables.  'base' option is slow but verbose.  'dplyr' is fast and quiet.
  + _find_optimal_components()_ Finds optium number of components for PLSR.  'pls' package consists of choosing the model with fewest components that is still less than one standard error away from the overall best model. 'first min' consists of choosing the first component that gives statistically (t-test) the same result as the following component. This method finds the first 'plateau' in the PRESS diminution.  'firstMin' finds the first component that gives statistically (t-test) the same result as the overall best model.
  

3) Example Rmarkdown vignettes illustrating the various PLSR model fitting examples

### Linked dataset citations, DOIs, and EcoSIS IDs/URLs: <br>
1) Leaf reflectance plant functional gradient IFGG/KIT <br>
Target variable: LMA <br>
EcoSIS URL: https://ecosis.org/package/leaf-reflectance-plant-functional-gradient-ifgg-kit <br>
EcoSIS ID: 3cf6b27e-d80e-4bc7-b214-c95506e46daa <br>
Rpubs example output: https://rpubs.com/sserbin/661964

2) Fresh Leaf Spectra to Estimate LMA over NEON domains in eastern United States <br>
Target variable: LMA <br>
EcoSIS URL: https://ecosis.org/package/fresh-leaf-spectra-to-estimate-lma-over-neon-domains-in-eastern-united-states <br>
EcoSIS ID: 5617da17-c925-49fb-b395-45a51291bd2d <br>
Rpubs example output: https://rpubs.com/sserbin/661976

3) Leaf spectra of 36 species growing in Rosa rugosa invaded coastal grassland communities in Belgium <br>
Target variable: LMA, leaf nitrogen <br>
EcoSIS URL: https://ecosis.org/package/leaf-spectra-of-36-species-growing-in-rosa-rugosa-invaded-coastal-grassland-communities-in-belgium <br>
EcoSIS ID: 9db4c5a2-7eac-4e1e-8859-009233648e89 <br>
DOI: https://doi.org/doi:10.21232/9nr6-sq54
Rpubs LMA example output: https://rpubs.com/sserbin/661963 <br>
Rpubs LeafN example output: https://rpubs.com/sserbin/661958

### Article citation:
TBD

### Script authors:
Julien Lamour, Jeremiah Anderson, Ken Davidson, Shawn P. Serbin [author list needs to be finalized]
