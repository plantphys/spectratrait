####################################################################################################
#
#  
#   An expanded example "How-to" script illustrating the use of PLSR modeling to develop a 
#   spectra-trait algorithm to estimate leaf mass area with leaf-level spectroscopy data. The 
#   example is built from published data source from the EcoSIS spectral database. This examples
#   illustrates an approach to quantify model prediction uncertainty based on a jackknife analysis
#
#   Spectra and trait data source:
#   https://ecosis.org/package/fresh-leaf-spectra-to-estimate-lma-over-neon-domains-in-eastern-united-states
#
#    Notes:
#    * The author notes the code is not the most elegant or clean, but is functional 
#    * Questions, comments, or concerns can be sent to sserbin@bnl.gov
#    * Code is provided under GNU General Public License v3.0 
#
####################################################################################################


#--------------------------------------------------------------------------------------------------#
### Install and load required R packages
list.of.packages <- c("devtools","readr","RCurl","httr","pls","dplyr","reshape2",
                      "ggplot2","gridExtra")  # packages needed for script
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load libraries
invisible(lapply(list.of.packages, library, character.only = TRUE))
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Setup other functions and options
github_dir <- file.path("~/Data/GitHub/PLSR_for_plant_trait_prediction/R_Scripts/")
source_from_gh <- FALSE
if (source_from_gh) {
  # Source helper functions from GitHub
  devtools::source_url("https://raw.githubusercontent.com/TESTgroup-BNL/PLSR_for_plant_trait_prediction/master/R_Scripts/functions.R")
} else {
  functions <- file.path(github_dir,"functions.R")
  source(functions)
}

# not in
`%notin%` <- Negate(`%in%`)

# Script options
pls.options(plsralg = "oscorespls")
pls.options("plsralg")
pls.options()$parallel

# Default par options
opar <- par(no.readonly = T)


# What is the target variable?
inVar <- "LMA_gDW_m2"

# What is the source dataset from EcoSIS?
ecosis_id <- "5617da17-c925-49fb-b395-45a51291bd2d"
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Set working directory (scratch space)
outdir <- tempdir()
setwd(outdir) # set working directory
getwd()  # check wd
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Get source dataset from EcoSIS
get_ecosis_data(ecosis_id = ecosis_id)
#--------------------------------------------------------------------------------------------------#




