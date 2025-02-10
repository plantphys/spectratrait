####################################################################################################
#
#  
#    Notes:
#    * The author notes the code is not the most elegant or clean, but is functional 
#    * Questions, comments, or concerns can be sent to sserbin@bnl.gov
#    * Code is provided under GNU General Public License v3.0 
#
####################################################################################################


#--------------------------------------------------------------------------------------------------#
### Load libraries
list.of.packages <- c("pls","dplyr","reshape2","here","plotrix","ggplot2","gridExtra",
                      "spectratrait")
invisible(lapply(list.of.packages, library, character.only = TRUE))
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Setup options

# Default par options
opar <- par(no.readonly = T)

# Specify output directory, output_dir 
# Options: 
# tempdir - use a OS-specified temporary directory 
# user defined PATH - e.g. "~/scratch/PLSR"
output_dir <- "tempdir"
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Set ecosis dataset
# https://ecosis.org/package/fresh-leaf-spectra-to-estimate-leaf-traits-for-california-ecosystems
# title: Fresh Leaf Spectra to Estimate Leaf Traits for California Ecosystems

# What is the target variable?
inVar <- "Nmass_g_g" # unclear if N in this dataset is Nmass or Narea. Assuming Nmass

# What is the source dataset from EcoSIS?
ecosis_id <- "0fadcc45-f79e-4fd3-a6ca-8afaf26ae299"
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
### Grab ecosis data
print(paste0("Output directory: ",getwd()))  # check wd
dat_raw <- spectratrait::get_ecosis_data(ecosis_id = ecosis_id)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### PLSR Coefficients - Grab from GitHub
git_repo <- "https://raw.githubusercontent.com/serbinsh/NASA_FFT_Leaf_Spectra-Trait_Models/refs/heads/master/"
print("**** Downloading PLSR coefficients ****")
githubURL <- paste0(git_repo,"PLSR_model_coefficients/leaf_Nitrogen/FFT_Leaf_Nitrogen_PLSR_Coefficients_11comp.csv")
LeafN.plsr.coeffs <- spectratrait::source_GitHubData(githubURL)
rm(githubURL)
githubURL <- paste0(git_repo,"PLSR_model_coefficients/leaf_Nitrogen/FFT_Leaf_Nitrogen_Jackkife_PLSR_Coefficients.csv")
LeafN.plsr.jk.coeffs <- spectratrait::source_GitHubData(githubURL)
rm(githubURL)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Setup target spectral and trait data
### Create plsr dataset
Start.wave <- 1500
End.wave <- 2400
wv <- seq(Start.wave,End.wave,1)
Spectra <- as.matrix(dat_raw[,names(dat_raw) %in% wv])
colnames(Spectra) <- c(paste0("Wave_",wv))
sample_info <- dat_raw[,names(dat_raw) %notin% seq(350,2500,1)]
head(sample_info)

sample_info2 <- sample_info %>%
  select(Plant_Species=`species`, Nmass_g_g=`Leaf nitrogen content per leaf dry mass`)
head(sample_info2)

plsr_data <- data.frame(sample_info2,Spectra)
rm(sample_info,sample_info2,Spectra)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
#### Example data cleaning.  End user needs to do what's appropriate for their 
#### data.  This may be an iterative process.
# Keep only complete rows of inVar and spec data before fitting
plsr_data <- plsr_data[complete.cases(plsr_data[,names(plsr_data) %in% 
                                                  c(inVar,paste0("Wave_",wv))]),]
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
print("**** Applying PLSR model to estimate leaf N mass from spectral observations ****")
# setup model
dims <- dim(LeafN.plsr.coeffs)
LeafN.plsr.intercept <- LeafN.plsr.coeffs[1,]
LeafN.plsr.coeffs <- data.frame(LeafN.plsr.coeffs[2:dims[1],])
names(LeafN.plsr.coeffs) <- c("wavelength","coefs")
LeafN.plsr.coeffs.vec <- as.vector(LeafN.plsr.coeffs[,2])
sub_spec <- droplevels(plsr_data[,which(names(plsr_data) %in% 
                                          paste0("Wave_",seq(Start.wave,End.wave,1)))])
#sub_spec <- sub_spec*0.01 # convert to 0-1
plsr_pred <- as.matrix(sub_spec) %*% LeafN.plsr.coeffs.vec + LeafN.plsr.intercept[,2]
leafN <- plsr_pred[,1]
names(leafN) <- "PLSR_LeafN_g_g"

# organize output
LeafN.PLSR.dataset <- data.frame(plsr_data[,which(names(plsr_data) %notin% 
                                                      paste0("Wave_",seq(Start.wave,End.wave,1)))],
                                   PLSR_LeafN_g_g=leafN, PLSR_Residuals=leafN-plsr_data[,inVar])
head(LeafN.PLSR.dataset)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
print("**** Generate PLSR uncertainty estimates ****")
jk_coef <- data.frame(LeafN.plsr.jk.coeffs[,3:dim(LeafN.plsr.jk.coeffs)[2]])
jk_coef <- t(jk_coef)
head(jk_coef)[,1:6]
jk_int <- t(LeafN.plsr.jk.coeffs[,2])
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
sd_res <- sd(LeafN.PLSR.dataset$PLSR_Residuals)
sd_tot <- sqrt(sd_mean^2+sd_res^2)
LeafN.PLSR.dataset$LCI <- Interval_Conf[1,]
LeafN.PLSR.dataset$UCI <- Interval_Conf[2,]
LeafN.PLSR.dataset$LPI <- LeafN.PLSR.dataset$PLSR_LeafN_g_g-1.96*sd_tot
LeafN.PLSR.dataset$UPI <- LeafN.PLSR.dataset$PLSR_LeafN_g_g+1.96*sd_tot
head(LeafN.PLSR.dataset)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
rmsep_percrmsep <- spectratrait::percent_rmse(plsr_dataset = LeafN.PLSR.dataset, 
                                              inVar = inVar, 
                                              residuals = LeafN.PLSR.dataset$PLSR_Residuals, 
                                              range="full")
RMSEP <- rmsep_percrmsep$rmse
perc_RMSEP <- rmsep_percrmsep$perc_rmse
r2 <- round(summary(lm(LeafN.PLSR.dataset$PLSR_LeafN_g_g~
                         LeafN.PLSR.dataset[,inVar]))$adj.r.squared,2)
expr <- vector("expression", 3)
expr[[1]] <- bquote(R^2==.(r2))
expr[[2]] <- bquote(RMSEP==.(round(RMSEP,2)))
expr[[3]] <- bquote("%RMSEP"==.(round(perc_RMSEP,2)))
rng_vals <- c(min(LeafN.PLSR.dataset$LPI), max(LeafN.PLSR.dataset$UPI))
par(mfrow=c(1,1), mar=c(4.2,5.3,1,0.4), oma=c(0, 0.1, 0, 0.2))
plotrix::plotCI(LeafN.PLSR.dataset$PLSR_LeafN_g_g,LeafN.PLSR.dataset[,inVar], 
                li=LeafN.PLSR.dataset$LPI, ui=LeafN.PLSR.dataset$UPI, gap=0.009,sfrac=0.000, 
                lwd=1.6, xlim=c(rng_vals[1], rng_vals[2]), ylim=c(rng_vals[1], rng_vals[2]), 
                err="x", pch=21, col="black", pt.bg=scales::alpha("grey70",0.7), scol="grey80",
                cex=2, xlab=paste0("Predicted ", paste(inVar), " (units)"),
                ylab=paste0("Observed ", paste(inVar), " (units)"),
                cex.axis=1.5,cex.lab=1.8)
abline(0,1,lty=2,lw=2)
plotrix::plotCI(LeafN.PLSR.dataset$PLSR_LeafN_g_g,LeafN.PLSR.dataset[,inVar], 
                li=LeafN.PLSR.dataset$LCI, ui=LeafN.PLSR.dataset$UCI, gap=0.009,sfrac=0.004, 
                lwd=1.6, xlim=c(rng_vals[1], rng_vals[2]), ylim=c(rng_vals[1], rng_vals[2]), 
                err="x", pch=21, col="black", pt.bg=scales::alpha("grey70",0.7), scol="black",
                cex=2, xlab=paste0("Predicted ", paste(inVar), " (units)"),
                ylab=paste0("Observed ", paste(inVar), " (units)"),
                cex.axis=1.5,cex.lab=1.8, add=T)
legend("topleft", legend=expr, bty="n", cex=1.5)
legend("bottomright", legend=c("Prediction Interval","Confidence Interval"), 
       lty=c(1,1), col = c("grey80","black"), lwd=3, bty="n", cex=1.5)
box(lwd=2.2)
dev.copy(png,file.path(outdir,paste0(inVar,"_PLSR_Validation_Scatterplot.png")), 
         height=2800, width=3200,  res=340)
dev.off();
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
### EOF