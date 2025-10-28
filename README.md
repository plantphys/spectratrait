# PLSR modeling for the estimation of plant functional traits
This repository contains source code and example scripts to illustrate best-practices for fitting, evaluating, and reporting spectra-trait PLSR models. This code and example scripts encompass several possibilities that you may encounter when carrying out PLSR model fitting. Start by reading *Burnett et al. (2021)*, then work through the scripts or vignettes  we have provided here to get experience in developing models for estimating plant functional traits using spectral measurements!   You can also explore examples available in the [Wiki](https://github.com/plantphys/spectratrait/wiki) pages

### Article citation:
Burnett AC, Anderson J, Davidson KD, Ely KS, Lamour J, Li Q, Morrison BD, Yang D, Rogers A, Serbin SP (2021) A best-practice guide to predicting plant traits from leaf-level hyperspectral data using partial least squares regression. Journal of Experimental Botany. https://doi.org/10.1093/jxb/erab295

### Source code citation:
[![DOI](https://zenodo.org/badge/222699149.svg)](https://zenodo.org/doi/10.5281/zenodo.4330119)

### EcoSML
https://ecosml.org/package/github/TESTgroup-BNL/spectratrait

### Getting started, tips and tricks:
* If you are new to R you should start by reading https://support.rstudio.com/hc/en-us/articles/201141096-Getting-Started-with-R & https://www.dataquest.io/blog/tutorial-getting-started-with-r-and-rstudio/
* Software requirements: R software (version 4.0 or above) and preferred operating environment (e.g. RStudio). 
* Install package dependencies and the spectratrait package: See the Depends and INSTALL sections below
* To work with the repository locally, clone the repository to your local machine (https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository). Once you have the repository on your local machine you can run the scripts in inst/scripts or vignettes folders. You can also start editing the code yourself or contributing to the development of the package through new pull requests (https://guides.github.com/activities/hello-world/)
* Or if you don't want to obtain the code via cloning the repository, you can click the green "Code" button at the top of this page and select "Download ZIP". Extract the contents of the ZIP in your preferred location. Use RStudio to open your selected example script from the inst/scripts folder and then run or source the code.

### Depends: 
ggplot2 (>= 3.3.2), remotes (>= 2.2.0), devtools (>= 2.3.1), readr (>= 1.3.1), RCurl (>= 1.98-1.2), 
httr (>= 1.4.2), pls (>= 2.7-2), magrittr (>= 2.0.1), dplyr (>= 1.0.1), reshape2 (>= 1.4.4), here (>= 0.1), 
plotrix (>= 3.7-8), gridExtra (>= 2.3), scales (>= 1.1.1), knitr (>= 1.4.2)

### INSTALL
spectratrait is not currently on CRAN, but you can install from GitHub using devtools().  First, make sure you have all of the package dependencies installed. You can do this either by 1) installing the packages individually using install.packages(), for example:

``` r
install.packages("pls")
install.packages("ggplot2")
...
```

and so forth until all of the dependencies (listed above in the "Depends" section) are installed. **Note** - you should pay careful attention at this stage to any R messages in your terminal alerting you that you need to update existing or install new R packages. These messages usually show up after you attempt to run install.packages() and require you 
to respond in your terminal to a y/n or multiple choice question before the install can continue.

Or 2) you can also run or source the "install_dependencies.R" script located in inst/scripts which should also install all of the required dependencies.  **Note** - again you will need to watch for any R prompts to update packages in order for the install to proceed correctly.

Finally, to complete the installation you will also need to install the spectratrait package itself.  You can do this by copying and pasting the command below into your R or RStudio (preferred) terminal.

``` r
# to install the master branch version
devtools::install_github(repo = "plantphys/spectratrait", dependencies=TRUE)

# to install the master branch version - with Vignettes (though slower)
devtools::install_github(repo = "plantphys/spectratrait", dependencies=TRUE, build_vignettes = TRUE)

# to install a specific release, for example release 1.0.5
devtools::install_github(repo = "plantphys/spectratrait@v1.0.5", dependencies=TRUE)

# or a specific branch, e.g. a branch named devbranch
devtools::install_github(repo = "plantphys/spectratrait", ref = "devbranch", dependencies=TRUE)
```

## Contains:
1. Core package functions are located in the in the main "R" folder
2. inst/scripts contains example PLSR workflows for fitting example leaf and canopy spectra-trait PLSR models for different leaf traits, including LMA and foliar    nitrogen
3. Example datasets that can be loaded in your R environment using the base load() function can be found in the data/ folder
4. man - the manual pages that are accessible in R 
5. tests - package tests to check that functions are still operational and produce the expected results
6. vignettes - example Rmarkdown and github markdown vignettes illustrating the various PLSR model fitting examples. These can be used to learn how to use the      PLSR workflow and associated functions for new applications
7. spectratrait_X.X.X.pdf (where X.X.X is the current release number) is the pdf documentation

### Linked dataset citations, DOIs, and EcoSIS IDs/URLs: <br>
1) Leaf spectra, structural and biochemical leaf traits of eight crop species (Ely et al., 2019) <br>
EcoSIS URL: https://ecosis.org/package/leaf-spectra--structural-and-biochemical-leaf-traits-of-eight-crop-species <br>
EcoSIS ID: 25770ad9-d47c-428b-bf99-d1543a4b0ec9 <br>
DOI: https://doi.org/doi:10.21232/C2GM2Z <br>
Rpubs LeafN bootstrap example output: https://rpubs.com/sserbin/spectratrait_ex1 <br>
Rpubs LeafN bootstrap by group (species) example output: https://rpubs.com/sserbin/spectratrait_ex2 <br>

2) Leaf reflectance plant functional gradient IFGG/KIT <br>
Target variable: SLA <br>
EcoSIS URL: https://ecosis.org/package/leaf-reflectance-plant-functional-gradient-ifgg-kit <br>
EcoSIS ID: 3cf6b27e-d80e-4bc7-b214-c95506e46daa <br>
Rpubs example output: https://rpubs.com/sserbin/spectratrait_ex3 <br>

3) Fresh leaf spectra to estimate LMA over NEON domains in eastern United States <br>
Target variable: LMA <br>
EcoSIS URL: https://ecosis.org/package/fresh-leaf-spectra-to-estimate-lma-over-neon-domains-in-eastern-united-states <br>
EcoSIS ID: 5617da17-c925-49fb-b395-45a51291bd2d <br>
DOI: https://doi.org/doi:10.21232/9831-rq60 <br>
Rpubs example output: https://rpubs.com/sserbin/spectratrait_ex4 <br>
Rpubs example showing Serbin et al. (2019) applied to NEON data: https://rpubs.com/sserbin/spectratrait_ex9 <br>

4) Canopy spectra to map foliar functional traits over NEON domains in eastern United States <br>
Target variable: leaf nitrogen <br>
EcoSIS URL: https://ecosis.org/package/canopy-spectra-to-map-foliar-functional-traits-over-neon-domains-in-eastern-united-states <br>
EcoSIS ID: b9dbf3db-5b9c-4ab2-88c2-26c8b39d0903 <br>
DOI: https://doi.org/doi:10.21232/e2jt-5209 <br>
Rpubs leaf nitrogen example output: https://rpubs.com/sserbin/spectratrait_ex5 <br>

5) Leaf spectra of 36 species growing in Rosa rugosa invaded coastal grassland communities in Belgium <br>
Target variable: LMA, leaf nitrogen <br>
EcoSIS URL: https://ecosis.org/package/leaf-spectra-of-36-species-growing-in-rosa-rugosa-invaded-coastal-grassland-communities-in-belgium <br>
EcoSIS ID: 9db4c5a2-7eac-4e1e-8859-009233648e89 <br>
DOI: https://doi.org/doi:10.21232/9nr6-sq54 <br>
Rpubs LeafN example output: https://rpubs.com/sserbin/spectratrait_ex6 <br>
Rpubs LeafN bootstrap example output: https://rpubs.com/sserbin/spectratrait_ex7 <br>
Rpubs LMA example output: https://rpubs.com/sserbin/spectratrait_ex8 <br>

## Build status
Auto-run PLSR example:
[![.github/workflows/run_plsr_example_auto.yaml](https://github.com/plantphys/spectratrait/actions/workflows/run_plsr_example_auto.yaml/badge.svg?branch=main)](https://github.com/plantphys/spectratrait/actions/workflows/run_plsr_example_auto.yaml) <br>
CI run PLSR example:
[![ci-run_PLSR_example](https://github.com/plantphys/spectratrait/actions/workflows/ci-run_plsr_example.yaml/badge.svg?branch=main)](https://github.com/plantphys/spectratrait/actions/workflows/ci-run_plsr_example.yaml) <br>
CI OS and R Release Checks:
[![R-CMD-check-OS-R](https://github.com/plantphys/spectratrait/actions/workflows/check-os.yaml/badge.svg?branch=main)](https://github.com/plantphys/spectratrait/actions/workflows/check-os.yaml) <br>
Weekly CI Checks:
[![R-CMD-check-Weekly](https://github.com/plantphys/spectratrait/actions/workflows/ci-weekly.yaml/badge.svg?branch=main)](https://github.com/plantphys/spectratrait/actions/workflows/ci-weekly.yaml) <br>
EcoSIS API Check:
[![run_ecosis_pull_example](https://github.com/plantphys/spectratrait/actions/workflows/run_ecosis_pull_example.yaml/badge.svg?branch=main)](https://github.com/plantphys/spectratrait/actions/workflows/run_ecosis_pull_example.yaml)
