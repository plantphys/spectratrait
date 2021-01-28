####################################################################################################
#
#  
#    Notes:
#    * The author notes the code is not the most elegant or clean, but is functional 
#    * Questions, comments, or concerns can be sent to sserbin@bnl.gov
#    * Code is provided under GNU General Public License v3.0 
#
####################################################################################################

rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files

#--------------------------------------------------------------------------------------------------#
### Install and load required R packages
list.of.packages <- list.of.packages <- c("readr","scales","plotrix","httr","devtools",
                                          "dplyr","here")
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=c("Depends", "Imports",
                                                                       "LinkingTo"))
version_requirements <- c("3.3.2")
if (!packageVersion("ggplot2") >= version_requirements[1]) {
  remotes::install_version(package="ggplot2", version=paste0(">= ", version_requirements), 
                           dependencies=c("Depends", "Imports", "LinkingTo"), upgrade="ask",
                           quiet=TRUE)
}
# Load libraries
invisible(lapply(list.of.packages, library, character.only = TRUE))
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Setup other functions and options
github_dir <- file.path(here(),"R_Scripts")
source_from_gh <- TRUE
if (source_from_gh) {
  # Source helper functions from GitHub
  print("*** GitHub hash of functions.R file:")
  devtools::source_url("https://raw.githubusercontent.com/TESTgroup-BNL/PLSR_for_plant_trait_prediction/master/R_Scripts/functions.R")
} else {
  functions <- file.path(github_dir,"functions.R")
  source(functions)
}

# not in
`%notin%` <- Negate(`%in%`)

# Default par options
opar <- par(no.readonly = T)

# What is the target variable?
inVar <- "LMA_gDW_m2"

# What is the source dataset from EcoSIS?
ecosis_id <- "5617da17-c925-49fb-b395-45a51291bd2d"

# Specify output directory, output_dir 
# Options: 
# tempdir - use a OS-specified temporary directory 
# user defined PATH - e.g. "~/scratch/PLSR"
output_dir <- "tempdir"
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Set working directory
if (output_dir=="tempdir") {
  outdir <- tempdir()
} else {
  if (! file.exists(output_dir)) dir.create(output_dir,recursive=TRUE)
  outdir <- file.path(path.expand(output_dir))
}
setwd(outdir) # set working directory
getwd()  # check wd
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### PLSR Coefficients - Grab from GitHub
git_repo <- "https://raw.githubusercontent.com/serbinsh/SSerbin_etal_2019_NewPhytologist/master/"
print("**** Downloading PLSR coefficients ****")
githubURL <- paste0(git_repo,"SSerbin_multibiome_lma_plsr_model/sqrt_LMA_gDW_m2_PLSR_Coefficients_10comp.csv")
LeafLMA.plsr.coeffs <- source_GitHubData(githubURL)
rm(githubURL)
githubURL <- paste0(git_repo,"SSerbin_multibiome_lma_plsr_model/sqrt_LMA_gDW_m2_Jackkife_PLSR_Coefficients.csv")
LeafLMA.plsr.jk.coeffs <- source_GitHubData(githubURL)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Get source dataset from EcoSIS
dat_raw <- get_ecosis_data(ecosis_id = ecosis_id)
head(dat_raw)
names(dat_raw)[1:40]
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Prepare new data for estimation
Start.wave <- 500
End.wave <- 2400
wv <- seq(Start.wave,End.wave,1)
Spectra <- as.matrix(dat_raw[,names(dat_raw) %in% wv])
colnames(Spectra) <- c(paste0("Wave_",wv))
head(Spectra)[1:6,1:10]
sample_info <- dat_raw[,names(dat_raw) %notin% seq(350,2500,1)]
head(sample_info)

sample_info2 <- sample_info %>%
  select(Domain,Functional_type,Sample_ID,USDA_Species_Code=`USDA Symbol`,LMA_gDW_m2=LMA)
head(sample_info2)

plsr_data <- data.frame(sample_info2,Spectra)
rm(sample_info,sample_info2,Spectra)
dim(plsr_data)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
#### Example data cleaning.  End user needs to do what's appropriate for their 
#### data.  This may be an iterative process.
# Keep only complete rows of inVar and spec data before fitting
plsr_data <- plsr_data[complete.cases(plsr_data[,names(plsr_data) %in% 
                                                  c(inVar,paste0("Wave_",wv))]),]
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
print("**** Applying PLSR model to estimate LMA from spectral observations ****")
# setup model
dims <- dim(LeafLMA.plsr.coeffs)
LeafLMA.plsr.intercept <- LeafLMA.plsr.coeffs[1,]
LeafLMA.plsr.coeffs <- data.frame(LeafLMA.plsr.coeffs[2:dims[1],])
names(LeafLMA.plsr.coeffs) <- c("wavelength","coefs")
LeafLMA.plsr.coeffs.vec <- as.vector(LeafLMA.plsr.coeffs[,2])
sub_spec <- droplevels(plsr_data[,which(names(plsr_data) %in% 
                                                   paste0("Wave_",seq(Start.wave,End.wave,1)))])
plsr_pred <- as.matrix(sub_spec) %*% LeafLMA.plsr.coeffs.vec + LeafLMA.plsr.intercept[,2]
leafLMA <- plsr_pred[,1]^2  # convert to standard LMA units from sqrt(LMA)
names(leafLMA) <- "PLSR_LMA_gDW_m2"

# organize output
LeafLMA.PLSR.dataset <- data.frame(plsr_data[,which(names(plsr_data) %notin% 
                                                      paste0("Wave_",seq(Start.wave,End.wave,1)))],
                                   PLSR_LMA_gDW_m2=leafLMA, PLSR_Residuals=leafLMA-plsr_data[,inVar])
head(LeafLMA.PLSR.dataset)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
print("**** Generate PLSR uncertainty estimates ****")
jk_coef <- data.frame(LeafLMA.plsr.jk.coeffs[,3:dim(LeafLMA.plsr.jk.coeffs)[2]])
jk_coef <- t(jk_coef)
head(jk_coef)[,1:6]
jk_int <- t(LeafLMA.plsr.jk.coeffs[,2])
head(jk_int)[,1:6]

jk_pred <- as.matrix(sub_spec) %*% jk_coef + matrix(rep(jk_int, length(plsr_data[,inVar])), 
                                         byrow=TRUE, ncol=length(jk_int))
jk_pred <- jk_pred^2
head(jk_pred)[,1:6]
dim(jk_pred)
interval <- c(0.025,0.975)
Interval_Conf <- apply(X = jk_pred, MARGIN = 1, FUN = quantile, 
                       probs=c(interval[1], interval[2]))
sd_mean <- apply(X = jk_pred, MARGIN = 1, FUN =sd)
sd_res <- sd(LeafLMA.PLSR.dataset$PLSR_Residuals)
sd_tot <- sqrt(sd_mean^2+sd_res^2)
LeafLMA.PLSR.dataset$LCI <- Interval_Conf[1,]
LeafLMA.PLSR.dataset$UCI <- Interval_Conf[2,]
LeafLMA.PLSR.dataset$LPI <- LeafLMA.PLSR.dataset$PLSR_LMA_gDW_m2-1.96*sd_tot
LeafLMA.PLSR.dataset$UPI <- LeafLMA.PLSR.dataset$PLSR_LMA_gDW_m2+1.96*sd_tot
head(LeafLMA.PLSR.dataset)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
RMSEP <- sqrt(mean(LeafLMA.PLSR.dataset$PLSR_Residuals^2))
pecr_RMSEP <- RMSEP/mean(LeafLMA.PLSR.dataset[,inVar])*100
r2 <- round(summary(lm(LeafLMA.PLSR.dataset$PLSR_LMA_gDW_m2~LeafLMA.PLSR.dataset[,inVar]))$adj.r.squared,2)
expr <- vector("expression", 3)
expr[[1]] <- bquote(R^2==.(r2))
expr[[2]] <- bquote(RMSEP==.(round(RMSEP,2)))
expr[[3]] <- bquote("%RMSEP"==.(round(pecr_RMSEP,2)))
rng_vals <- c(min(LeafLMA.PLSR.dataset$LPI), max(LeafLMA.PLSR.dataset$UPI))
par(mfrow=c(1,1), mar=c(4.2,5.3,1,0.4), oma=c(0, 0.1, 0, 0.2))
plotCI(LeafLMA.PLSR.dataset$PLSR_LMA_gDW_m2,LeafLMA.PLSR.dataset[,inVar], 
       li=LeafLMA.PLSR.dataset$LCI, ui=LeafLMA.PLSR.dataset$UCI, gap=0.009,sfrac=0.004, 
       lwd=1.6, xlim=c(rng_vals[1], rng_vals[2]), ylim=c(rng_vals[1], rng_vals[2]), 
       err="x", pch=21, col="black", pt.bg=alpha("grey70",0.7), scol="grey80",
       cex=2, xlab=paste0("Predicted ", paste(inVar), " (units)"),
       ylab=paste0("Observed ", paste(inVar), " (units)"),
       cex.axis=1.5,cex.lab=1.8)
abline(0,1,lty=2,lw=2)
legend("topleft", legend=expr, bty="n", cex=1.5)
box(lwd=2.2)
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
### EOF