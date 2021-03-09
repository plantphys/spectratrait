# PLSR modeling for the estimation of plant functional traits
This repository contains example scripts illustrating best-practices for fitting, evaluating, and reporting leaf-level spectra-trait PLSR models. These scripts encompass several possibilities that you may encounter doing PLSR. Start by reading *Burnett et al. in review*, then work through the scripts or vignettes.

### Article citation:
Burnett AC, Anderson J, Davidson KD, Ely KS, Lamour J, Li Q, Morrison BD, Yang D, Rogers A, Serbin SP (in review) A best-practice guide to predicting plant traits from leaf-level hyperspectral data using partial least squares regression. Journal of Experimental Botany. In Review.

### Source code citation:
[![DOI](https://zenodo.org/badge/222699149.svg)](https://zenodo.org/badge/latestdoi/222699149)

### EcoSML
https://ecosml.org/package/github/TESTgroup-BNL/PLSR_for_plant_trait_prediction

### Getting started, tips and tricks:
* If you are new to R you should start by reading https://support.rstudio.com/hc/en-us/articles/201141096-Getting-Started-with-R & https://www.dataquest.io/blog/tutorial-getting-started-with-r-and-rstudio/
* Software requirements: R software (version 4.0 or above) and preferred operating environment (e.g. RStudio). 
* Install package dependencies and the spectratrait package: See the Depends and INSTALL sections below
* To work with the repository locally, clone the repository to your local machine (https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository). Once you have the repository on your local machine you can run the scripts in inst/scripts or vignettes folders. You can also start editing the code yourself or contributing to the development of the package through new pull requests (https://guides.github.com/activities/hello-world/)
* Or if you don't want to obtain the code via cloning the repository, you can click the green "Code" button at the top of this page and select "Download ZIP". Extract the contents of the ZIP in your preferred location. Use RStudio to open your selected example script from the inst/scripts folder and then run or source the code.

### Depends: 
ggplot2 (>= 3.3.2), remotes (>= 2.2.0), devtools (>= 2.3.1), readr (>= 1.3.1), RCurl (>= 1.98-1.2), 
httr (>= 1.4.2), pls (>= 2.7-2), dplyr (>= 1.0.1), reshape2 (>= 1.4.4), here (>= 0.1), 
plotrix (>= 3.7-8), gridExtra (>= 2.3), scales (>= 1.1.1)

### INSTALL
spectratrait is not currently on CRAN, but you can install from GitHub using devtools().  First, make sure you have all of the package dependencies installed. You can do this either by 1) installing the packages individually using install.packages(), for example:

``` r
install.packages("pls")
install.packages("ggplot2")
...
```
Or 2) you can install all the packages at once as a list, though you need to carefully pay attention to any R messages about needing to uppdate existing or install new R packages at the same time. This usually shows up as a message with one or more packages that R would like to update requiring you to respond before install can continue

``` r
install.packages(c("devtools","remotes","readr","RCurl","httr","pls","dplyr","reshape2",
"here","plotrix","scales","ggplot2","gridExtra"))
```

You can copy and paste the above line directly into your R or RStudio terminal. Finally, you can also run or source the "install_dependencies.R" script located in inst/scripts which should also install the dependencies.  Note that again you will need to watch for any R prompts to update pacakges in order for install to proceed correctly.

Finally, to complete install you will also need to install the spectratrait package itself.  You can do this by copying and pasting the command below into your R or RStudio (preferred) terminal.

``` r
# to install the master branch version
devtools::install_github(repo = "TESTgroup-BNL/PLSR_for_plant_trait_prediction", 
dependencies=TRUE)

# to install a specific release, for example release 1.0.0
devtools::install_github(repo = "TESTgroup-BNL/PLSR_for_plant_trait_prediction@v1.0.0", 
dependencies=TRUE)

# or a specific branch, e.g. a branch named devbranch
devtools::install_github(repo = "TESTgroup-BNL/PLSR_for_plant_trait_prediction", 
ref = "devbranch", dependencies=TRUE)
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
1) Leaf reflectance plant functional gradient IFGG/KIT <br>
Target variable: SLA <br>
EcoSIS URL: https://ecosis.org/package/leaf-reflectance-plant-functional-gradient-ifgg-kit <br>
EcoSIS ID: 3cf6b27e-d80e-4bc7-b214-c95506e46daa <br>
Rpubs example output: https://rpubs.com/sserbin/722040

2) Fresh Leaf Spectra to Estimate LMA over NEON domains in eastern United States <br>
Target variable: LMA <br>
EcoSIS URL: https://ecosis.org/package/fresh-leaf-spectra-to-estimate-lma-over-neon-domains-in-eastern-united-states <br>
EcoSIS ID: 5617da17-c925-49fb-b395-45a51291bd2d <br>
DOI: https://doi.org/doi:10.21232/9831-rq60 <br>
Rpubs example output: https://rpubs.com/sserbin/722017

3) Leaf spectra of 36 species growing in Rosa rugosa invaded coastal grassland communities in Belgium <br>
Target variable: LMA, leaf nitrogen <br>
EcoSIS URL: https://ecosis.org/package/leaf-spectra-of-36-species-growing-in-rosa-rugosa-invaded-coastal-grassland-communities-in-belgium <br>
EcoSIS ID: 9db4c5a2-7eac-4e1e-8859-009233648e89 <br>
DOI: https://doi.org/doi:10.21232/9nr6-sq54 <br>
Rpubs LMA example output: https://rpubs.com/sserbin/721898 <br>
Rpubs LeafN example output: https://rpubs.com/sserbin/721878 <br>
Rpubs LeafN bootstrap example output: https://rpubs.com/sserbin/721908

4) Leaf spectra, structural and biochemical leaf traits of eight crop species <br>
EcoSIS URL: https://ecosis.org/package/leaf-spectra--structural-and-biochemical-leaf-traits-of-eight-crop-species <br>
EcoSIS ID: 25770ad9-d47c-428b-bf99-d1543a4b0ec9 <br>
DOI: https://doi.org/doi:10.21232/C2GM2Z <br>

5) Canopy Spectra to Map Foliar Functional Traits over NEON domains in eastern United States <br>
Target variable: leaf nitrogen <br>
EcoSIS URL: https://ecosis.org/package/canopy-spectra-to-map-foliar-functional-traits-over-neon-domains-in-eastern-united-states <br>
EcoSIS ID: b9dbf3db-5b9c-4ab2-88c2-26c8b39d0903 <br>
DOI: https://doi.org/doi:10.21232/e2jt-5209 <br>
Rpubs leaf nitrogen example output: https://rpubs.com/sserbin/728124
