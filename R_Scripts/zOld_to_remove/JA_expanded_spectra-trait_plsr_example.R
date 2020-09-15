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
#
#    --- Last updated:  07.29.2020 By Shawn P. Serbin <sserbin@bnl.gov>
####################################################################################################


#--------------------------------------------------------------------------------------------------#
### Install and load required R packages
list.of.packages <- c("devtools", "pls", "parallel")  # packages needed for script
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load libraries
library(pls)

#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Setup other functions and options
# Source helper functions from GitHub
devtools::source_url("https://raw.githubusercontent.com/TESTgroup-BNL/How_to_PLSR/master/R_Scripts/functions.R")

# not in
`%notin%` <- Negate(`%in%`)

# Script options
pls.options(plsralg = "oscorespls")
pls.options("plsralg")
pls.options()$parallel
# NULL

# Default par options
opar <- par(no.readonly = T)

# What is the target variable?
inVar <- "LMA_gDW_m2"
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Set working directory (scratch space)
outdir <- tempdir()
setwd(outdir) # set working directory
getwd()  # check wd
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Get source dataset from EcoSIS

print("**** Downloading Ecosis data ****")
ecosis_id <- "5617da17-c925-49fb-b395-45a51291bd2d"  # NEON dataset
ecosis_file <- sprintf(
  "http://ecosis.org/api/package/%s/export?metadata=true",
  ecosis_id
)
message("Downloading data...")
dat_raw <- read.csv(ecosis_file, na.strings = "")
message("Download complete!")
head(dat_raw[,1:40])
names(dat_raw)[1:40]
rm(ecosis_file, ecosis_id)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Create full plsr dataset
Start_wave <- 500
End_wave <- 2400
wv <- seq(Start_wave,End_wave,1)
wv <- c(paste0("X", wv))
spectra <- data.frame(dat_raw[,names(dat_raw) %in% wv])
names(spectra) <- gsub("X","Wave_",names(spectra))
head(spectra)[1:6,1:10]
sample_info <- dat_raw[ , -grep("X", names(dat_raw))]
head(sample_info)

sample_info2 <- sample_info[,c("Domain", "Functional_type", "Sample_ID", "USDA.Symbol", "LMA")]
names(sample_info2)[names(sample_info2) %in% "USDA.Symbol"] <- "USDA_Species_Code"
names(sample_info2)[names(sample_info2) %in% "LMA"] <- "LMA_gDW_m2"
head(sample_info2)

plsr_data <- data.frame(sample_info2,spectra)
head(plsr_data)[,1:10]
rm(sample_info,sample_info2,spectra, wv)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Create cal/val datasets
set.seed(2356812)

plsr_data$CalVal <- NA
split_var <- c("Domain", "USDA_Species_Code")
plsr_data$ID <- apply(plsr_data[, split_var], MARGIN = 1, FUN = function(x) paste(x, collapse = " "))
split_var_list <- unique(plsr_data$ID)
split_var_list <- split_var_list[ -grep("NA", split_var_list)]
prop <- 0.8

for(i in 1:length(split_var_list)){
  temp <- row.names(plsr_data[ plsr_data$ID == split_var_list[i], ])
  ## there should probably be more than 4 obs I'm guessing, so this may need adjusting
  if(length(temp) > 3){
    Cal <- sample(temp,round(prop*length(temp)))
    Val <- temp[!temp %in% Cal]
    plsr_data$CalVal[ row.names(plsr_data) %in% Cal ] <- "Cal"
    plsr_data$CalVal[ row.names(plsr_data) %in% Val ] <- "Val"
    p_cal <- length(Cal)/length(temp) * 100
    message(paste0(split_var_list[i], "   ", "Cal", ": ", p_cal, "%"))
  } else {
    message(paste(split_var_list[i], "Not enough observations"))
  }
}

plsr_data$ID <- NULL
# drop NA's  Not all observations had a species associated with it
plsr_data <- plsr_data[!is.na(plsr_data$CalVal), ]

rm(Cal, Val, split_var, split_var_list, p_cal, i, prop, temp)

print(paste("Cal observations: ",length(plsr_data$CalVal[plsr_data$CalVal == "Cal"])))
print(paste("Val observations: ",length(plsr_data$CalVal[plsr_data$CalVal == "Val"])))

par(mfrow=c(1,2), mar=c(4,4,3,0)+0.1)
nBins <- diff(range(plsr_data[,inVar]))/10
hist(plsr_data[plsr_data$CalVal == "Cal",inVar], breaks=nBins, main="Calibration", 
     xlab=inVar, ylab="Count")
hist(plsr_data[plsr_data$CalVal == "Val",inVar], breaks=nBins, main="Validation", 
     xlab=inVar, ylab="Count")
par(opar)

## write plsr data to CSV
write.csv(plsr_data,file=file.path(outdir,paste0(inVar,'_Full_PLSR_Dataset.csv')),row.names=FALSE)

#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Create calibration and validation PLSR datasets

spec <- as.matrix(plsr_data[,grep("Wave_", names(plsr_data)) ] )
plsr_data <- data.frame(plsr_data[,-grep("Wave", names(plsr_data))],Spectra=I(spec))
rm(spec)

# plot cal and val spectra
# summary stats on plots cut off, too far to the left.
par(mfrow=c(1,2)) # B, L, T, R
f.plot.spec(Z=plsr_data[plsr_data$CalVal=="Cal", "Spectra"],wv=seq(Start_wave,End_wave,1),plot_label="Calibration")
f.plot.spec(Z=plsr_data[plsr_data$CalVal=="Val", "Spectra"],wv=seq(Start_wave,End_wave,1),plot_label="Validation")
par(opar)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Use Jackknife permutation to determine optimal number of components
if(grepl("Windows", sessionInfo()$running)){
  pls.options(parallel =NULL)
} else {
  pls.options(parallel = parallel::detectCores()-1)
}

cal.plsr.data <- plsr_data[ plsr_data$CalVal == "Cal",]
nComps <- 20 # Not enough nComps sometimes, so changed to 20
iterations <- 20
seg <- 15
prop <- 0.70
jk.out <- matrix(data=NA,nrow=iterations,ncol=nComps) 
print("*** Running jacknife permutation test.  Please hang tight, this can take awhile ***")
start.time <- Sys.time()
for (i in 1:iterations) {
  message(paste("Running interation", i))
  rows <- sample(1:nrow(cal.plsr.data),floor(prop*nrow(cal.plsr.data)))
  sub.data <- cal.plsr.data[rows,]
  plsr.out <- plsr(as.formula(paste(inVar,"~","Spectra")), scale=FALSE, center=TRUE, ncomp=nComps, 
                   validation="CV", segments = seg, segment.type="interleaved", trace=FALSE, data=sub.data)
  resPRESS <- as.vector(plsr.out$validation$PRESS)
  jk.out[i,seq(plsr.out$validation$ncomp)]=resPRESS
}
end.time <- Sys.time()
end.time - start.time

boxplot(jk.out, xlab="Number of Components", ylab="PRESS")

results <- NULL
for(i in 1:(nComps-1)){
  p_value <- t.test(jk.out[,i], jk.out[,(i+1)])$p.value
  temp_results <- data.frame(Component=(i+1), P.value= round(p_value, 4))
  results <- rbind(results, temp_results)
}
results

# Simple final model validated with cross-validation.  Segmented cross-validation used
# given the very large sample size.  For models with fewer observations (e.g. <100) 
# LOO or leave-one-out cross validation is recommended

nComps <- min(results[results$P.value > 0.05, "Component"])
print(paste0("*** Optimal number of components based on t.test: ", nComps))

segs <- 30
#pls.options(parallel = NULL)
plsr.out <- plsr(as.formula(paste(inVar,"~","Spectra")),scale=FALSE,ncomp=nComps,validation="CV",
                 segments=segs, segment.type="interleaved",trace=TRUE,data=cal.plsr.data)
fit <- plsr.out$fitted.values[,1,nComps]
pls.options(parallel = NULL)

## Cleanup
rm(iterations, end.time, start.time, rows, seg, segs, temp_results, sub.data, jk.out,
   results, prop, resPRESS)

#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Generate some initial PLSR results
# External validation
val.plsr.data <- plsr_data[plsr_data$CalVal == "Val",]
par(mfrow=c(1,2)) # B, L, T, R
RMSEP(plsr.out, newdata = val.plsr.data)
plot(RMSEP(plsr.out,estimate=c("test"),newdata = val.plsr.data), main="MODEL RMSEP",
     xlab="Number of Components",ylab="Model Validation RMSEP",lty=1,col="black",cex=1.5,lwd=2)
box(lwd=2.2)

R2(plsr.out, newdata = val.plsr.data)
plot(R2(plsr.out,estimate=c("test"),newdata = val.plsr.data), main="MODEL R2",
     xlab="Number of Components",ylab="Model Validation R2",lty=1,col="black",cex=1.5,lwd=2)
box(lwd=2.2)
par(opar)

#calibration
cal.plsr.data$Fitted <- fit
cal.plsr.data$Residuals <- cal.plsr.data$Fitted - cal.plsr.data$LMA_gDW_m2

# validation
val.plsr.data$Fitted <- predict(plsr.out, newdata = val.plsr.data, ncomp=nComps, 
                                                 type="response")[,,1]
val.plsr.data$Residuals <- val.plsr.data$Fitted - val.plsr.data$LMA_gDW_m2

# Scatter setup
par(mfrow=c(2,2), mar=c(4,4,3,0)+0.1 )
plot_min <- min(cal.plsr.data$LMA_gDW_m2, cal.plsr.data$Fitted,
                val.plsr.data$LMA_gDW_m2, val.plsr.data$Fitted, na.rm = T)
plot_max <- max(cal.plsr.data$LMA_gDW_m2, cal.plsr.data$Fitted,
                val.plsr.data$LMA_gDW_m2, val.plsr.data$Fitted, na.rm = T)
plot_range <- c(plot_min, plot_max)
x_lab <- expression(paste("Predicted LMA (",g~m^{-2},")"))
y_lab <- expression(paste("Observed LMA (",g~m^{-2},")"))

# Cal Scatter
plsr_stats <- c(round(pls::R2(plsr.out)[[1]][nComps], 2),
                round(pls::RMSEP(plsr.out)[[1]][nComps], 2))

plot(plot_range, plot_range, xlab = x_lab, ylab = y_lab, type = "n", main="Calibration")
grid(NULL, NULL, lwd=1, lty=1)
points(cal.plsr.data$Fitted, cal.plsr.data$LMA_gDW_m2, pch=16)
abline(0, 1, lty=2, col = "grey", lwd=4)
mtext(bquote(R^2 == .(plsr_stats[1]) ~~ RMSE == .(plsr_stats[2])), line=0)

# Val Scatter
plsr_stats <- c(round(pls::R2(plsr.out, newdata = val.plsr.data)[[1]][nComps],2),
                round(pls::RMSEP(plsr.out, newdata = val.plsr.data)[[1]][nComps],2))

plot(plot_range, plot_range, xlab = x_lab, ylab = y_lab, type = "n", main="Validation")
grid(NULL, NULL, lwd=1, lty=1)
points(val.plsr.data$Fitted, val.plsr.data$LMA_gDW_m2, pch=16)
abline(0, 1, lty=2, col = "grey", lwd=4)
mtext(bquote(R^2 == .(plsr_stats[1]) ~~ RMSE == .(plsr_stats[2])), line=0)

# Cal Residuals Histogram
nBins <- round(diff(range(cal.plsr.data$Residuals))/.5, -1)
hist(cal.plsr.data$Residuals, breaks = nBins, border="lightgrey", xlab="Cal Residuals", ylab="Count", main=NULL )
abline(v=0, lty=2)

#Val Residuals Histogram
nBins <- round(diff(range(val.plsr.data$Residuals))/.5, -1)
hist(val.plsr.data$Residuals, breaks = nBins, border="lightgrey", xlab="Val Residuals", ylab="Count", main=NULL )
abline(v=0, lty=2)
par(opar)

#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### results by functional type and domain
# validation by functional type
plot_var <- "Functional_type"
plot_list <- unique(val.plsr.data[, plot_var])
col_list <- hcl.colors(length(plot_list), alpha=.75) ## can choose a different palette
#col_list <- rainbow(length(plot_list), alpha=0.75) ## also looks ok
plot(plot_range, plot_range, xlab = x_lab, ylab = y_lab, type = "n", main="Validation")
for(i in 1:length(plot_list)){
  plot_temp <- val.plsr.data[val.plsr.data[,plot_var] == plot_list[i], ]
  points(plot_temp$Fitted, plot_temp$LMA_gDW_m2, col = col_list[i], pch=16)
}
abline(0, 1, lty=2, col = "grey", lwd=4)
legend("topleft", legend=plot_list, col=col_list, bty="n", pch=16, horiz=T)

for(i in 1:length(plot_list)){
  plot_temp <- val.plsr.data[ val.plsr.data[plot_var] == plot_list[i], "Residuals"]
  hist_temp <- hist(plot_temp, round(diff(range(plot_temp))/.5, -1), plot=F)
  plot(hist_temp, col=col_list[i], xlim=range(val.plsr.data$Residuals),
       add=ifelse(i>1, T, F), border=col_list[i], ylab="Count", xlab="Residuals", main=NULL)
}
legend("topleft", legend=plot_list, col=col_list, bty="n", pch=16, horiz=T)

## By Domain
plot_var <- "Domain"
plot_list <- unique(val.plsr.data[, plot_var])
#col_list <- hcl.colors(length(plot_list), alpha=.75)
col_list <- rainbow(length(plot_list), alpha=0.75) ## works better for this one
plot(plot_range, plot_range, xlab = x_lab, ylab = y_lab, type = "n", main="Validation")
for(i in 1:length(plot_list)){
  plot_temp <- val.plsr.data[val.plsr.data[,plot_var] == plot_list[i], ]
  points(plot_temp$Fitted, plot_temp$LMA_gDW_m2, col = col_list[i], pch=16)
}
abline(0, 1, lty=2, col = "grey", lwd=4)
legend("topleft", legend=plot_list, col=col_list, bty="n", pch=16, horiz=T)

for(i in 1:length(plot_list)){
  plot_temp <- val.plsr.data[ val.plsr.data[plot_var] == plot_list[i], "Residuals"]
  hist_temp <- hist(plot_temp, round(diff(range(plot_temp))/.5, -1), plot=F)
  plot(hist_temp, col=col_list[i], xlim=range(val.plsr.data$Residuals),
       add=ifelse(i>1, T, F), border=col_list[i], ylab="Count", xlab="Residuals", main=NULL)
}
legend("topleft", legend=plot_list, col=col_list, bty="n", pch=16, horiz=T)

## Cleanup
rm(col_list, i, p_value, plot_list, plot_max, plot_min, plot_range, plot_temp, 
   plot_var, plsr_stats, x_lab, y_lab, hist_temp, nBins, fit)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Generate some useful outputs
cal.plsr.pred <- as.vector(plsr.out$fitted.values[,,nComps]) # Model fitted values. Predicted values
cal.plsr.CVpred <- as.vector(plsr.out$validation$pred[,,nComps]) # CV pred values
cal.CVresiduals <- as.vector(plsr.out$residuals[,,nComps]) # CV pred residuals
cal.output <- data.frame(cal.plsr.data[,which(names(cal.plsr.data) %notin% "Spectra")],
                         PLSR_Predicted=cal.plsr.pred, PLSR_CV_Predicted=cal.plsr.CVpred,
                         PLSR_CV_Residuals=cal.CVresiduals)
head(cal.output)
rm(cal.plsr.pred,cal.plsr.CVpred,cal.CVresiduals) #cleanup

predicted_val <- as.vector(predict(plsr.out, newdata = val.plsr.data, ncomp=nComps, type="response")[,,1])
predicted_val_residuals <- predicted_val-val.plsr.data[,inVar]
val.output <- data.frame(val.plsr.data[,which(names(cal.plsr.data) %notin% "Spectra")],
                         PLSR_Predicted=predicted_val,PLSR_Residuals=predicted_val_residuals)
head(val.output)
rm(predicted_val,predicted_val_residuals) #cleanup

coefs <- coef(plsr.out,ncomp=nComps,intercept=FALSE)
vips <- VIP(plsr.out)[nComps,]

# Coefficient and VIP plot for PLSR model !! This plotting could be improved !!
par(mfrow=c(2,1))
waves_list <- seq(Start_wave,End_wave,1)
plot(waves_list,coefs,cex=0.01,xlab="Wavelength (nm)",ylab="REG COEF")
lines(waves_list,coefs,lwd=2.5)
abline(h=0,lty=2, col = "grey", lwd=2)

plot(waves_list,vips,xlab="Wavelength (nm)",ylab="VIP",cex=0.01)
lines(waves_list,vips,lwd=3)
abline(h=0.8,lty=2, col = "grey", lwd=2)
par(opar)

#--------------------------------------------------------------------------------------------------#


#---------------- Export Model Output -------------------------------------------------------------#
print(paste("Output directory: ", getwd()))

# Observed versus predicted
write.csv(cal.output,file=file.path(outdir,paste0(inVar,'_Observed_PLSR_CV_Pred_',nComps,
                                                  'comp.csv')),row.names=FALSE)

# Validation data
write.csv(val.output,file=file.path(outdir,paste0(inVar,'_Val_PLSR_Pred_',nComps,
                                                  'comp.csv')),row.names=FALSE)

# Model coefficients
coefs_w_intercept <- coef(plsr.out,ncomp=nComps,intercept=TRUE)
write.csv(coefs_w_intercept,file=file.path(outdir,paste0(inVar,'_PLSR_Coefficients_',nComps,'comp.csv')),
          row.names=TRUE)

# PLSR VIP
write.csv(vips,file=file.path(outdir,paste0(inVar,'_PLSR_VIPs_',nComps,'comp.csv')))

# confirm files were written to temp space
print("**** PLSR output files: ")
list.files(getwd())[grep(pattern = inVar, list.files(getwd()))]

#--------------------------------------------------------------------------------------------------#


#---------------- Jackknife model evaluation ------------------------------------------------------#
nResamples_seq <- 1:100

pls.options(parallel =NULL)  ## actually slower using parallel on my setup (2.3 quad i7)

model_eval <- function(nResamples_seq, x, nComps){
  message(paste("Running iteration:", nResamples_seq))
  rows <- sample(1:nrow(x),floor(0.7*nrow(x)))
  cal.data.jk <- x[rows,]
  val.data.jk <- x[-rows,]

  ### Build PLSR model with training data
  # what is the scaled one for?  I don't see it used anywhere
  pls.jack <- plsr(as.formula(paste0(inVar,"~","Spectra")), scale=FALSE, ncomp=nComps, validation="none", data=cal.data.jk)
  pls.jack.scaled <- plsr(as.formula(paste0(inVar,"~","Spectra")), scale=TRUE, ncomp=nComps, validation="none", data=cal.data.jk)
  
  ### Estimate independent (Validation) samples
  plsr.val <- val.data.jk[,inVar]
  pred.val.data <- as.vector(predict(pls.jack,newdata=val.data.jk$Spectra,
                                     ncomp=nComps,type="response")[,,1])

  ### Coefficients and VIPs
  resam_coefs <- coef(pls.jack,ncomp=nComps,intercept=TRUE)
  coefs_scaled <- coef(pls.jack.scaled,ncomp=nComps,intercept=TRUE)
  vips <- c(NA, VIP(pls.jack)[nComps,])
  coef_vips <- as.data.frame(rbind(coefs=resam_coefs, coefs_scaled, vips))
  names(coef_vips)[1] <- "Intercept"
  coef_vips <- cbind(Interation = nResamples_seq, ID = row.names(coef_vips), coef_vips)
  
  # Error statistics
  Rsq.val <- R2(pls.jack,newdata = val.data.jk)$val[,,nComps+1] 
  val.residuals <- pred.val.data-plsr.val ## single value
  Val.bias <- mean(pred.val.data)-mean(plsr.val)
  MSEP <- mean(val.residuals^2)
  RMSEP.val <- sqrt(MSEP)
  PERC_RMSEP <- (RMSEP.val/(max(plsr.val)-min(plsr.val)))*100
  
  ### Store results of iteration
  Error_stats <- data.frame(Interation = nResamples_seq, ID = "Stats", Rsq=Rsq.val,
                            RMSEP=RMSEP.val, PERC_RMSEP=PERC_RMSEP, Bias=Val.bias)
  out <- merge(Error_stats, coef_vips, all=T)
 
  return(out)
}

eval_df <- do.call(rbind, lapply(nResamples_seq, model_eval, x=cal.plsr.data, nComps=nComps))


#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Histogram statistics
output_stats <- eval_df[eval_df$ID == "Stats", c("Rsq", "RMSEP", "PERC_RMSEP", "Bias")]
head(output_stats)

par(mfrow=c(1,3))
hist(output_stats$Rsq, breaks="FD", main="Resampled Rsq", ylab="Count", xlab=inVar)
hist(output_stats$RMSEP, breaks="FD", main="Resampled RMSEP", ylab="Count", xlab=inVar)
hist(output_stats$Bias, breaks="FD", main="Resampled Bias", ylab="Count", xlab=inVar)
par(opar)

#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### plot jackknife coefficients
plot_coefs <- eval_df[eval_df$ID == "coefs", grep("Wave_", names(eval_df))]
intercepts <-eval_df[eval_df$ID == "coefs", "Intercept"]
y_range <- range(plot_coefs, na.rm = T)

sum_stats_fun <- function(x){
  c(Mean = mean(x), SD=sd(x), Min=min(x), Max=max(x), Quant=quantile(x, c(0.025, 0.975)),
     Ttest_Pvalue = t.test(x)$p.value)
}
sum_stats <- sapply(cbind(intercepts, plot_coefs), sum_stats_fun)

## plot it
plot(range(waves_list), y_range, type="n", ylab="Coef", xlab="Wavelength (nm)")

polygon(c(waves_list ,rev(waves_list)),c(sum_stats[row.names(sum_stats)=="Max",-1], 
          rev(sum_stats[row.names(sum_stats)=="Min",-1])),col="grey50",border=NA)

polygon(c(waves_list ,rev(waves_list)), c(sum_stats[row.names(sum_stats)=="Quant.97.5%",-1], 
          rev(sum_stats[row.names(sum_stats)=="Quant.2.5%",-1])), col="grey70",border=NA)

lines(waves_list,coefs)

abline(h=0,lty=2,col="grey",lwd=1.5)
legend("bottomright",legend=c("Mean","95% CI", "Max/Min"),lty=c(1,1),
       col=c("black","grey70", "grey50"),lwd=3)

# P-vals
plot(waves_list, sum_stats[row.names(sum_stats) == "Ttest_Pvalue", -1], pch=21,
     bg="grey80", ylab="P-value", xlab="Wavelength (nm)", cex=2)

#--------------------------------------------------------------------------------------------------#


#---------------- Output jackknife results --------------------------------------------------------#
# JK Coefficents
write.csv(eval_df, file=file.path(outdir,paste0(inVar,'_Resample_PLSR_Coefficients_VIPS_Stats.csv')),
          row.names=FALSE)

write.csv(sum_stats, file=file.path(outdir,paste0(inVar,'_Jackkife_PLSR_Coefficient_Summary_Stats.csv')),
          row.names=TRUE)

#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# JK Val plot

coef_test_output <- NULL
intercept_coef <- eval_df[eval_df$ID=="coefs", c(grep("Wave_|Intercept", names(eval_df)))]
for(i in 1:nrow(intercept_coef)){
  temp <- val.plsr.data$Spectra %*% as.numeric(intercept_coef[i,-1])
  vals <- data.frame(rowSums(temp))+intercept_coef[i,1]
  coef_test_output <- cbind(coef_test_output, vals[,1])
}

pred_quant <- apply(coef_test_output,1,quantile,probs=c(0.025,0.975))

jk_val_plot_data <- data.frame(val.output, LL=pred_quant[1,], UL=pred_quant[2,])

plot_range <- range(jk_val_plot_data$PLSR_Predicted, jk_val_plot_data[,inVar])

plot(plot_range, plot_range, type="n", xlab=expression(Predicted ~ LMA ~ "(" ~g~m^-2~")"),
     ylab=expression(Observed ~ LMA ~ "(" ~g~m^-2~")"))
points(jk_val_plot_data$PLSR_Predicted, jk_val_plot_data[,inVar], pch=16)
arrows(x0=jk_val_plot_data$LL, x1=jk_val_plot_data$UL, y0=jk_val_plot_data$LMA_gDW_m2,
       length=0.05, angle=90, code=3)
abline(0,1, col="darkgrey", lty=2, lwd=3)

#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Output JK Coefficient test results
coef_test_output2 <- data.frame(coef_test_output)
names(coef_test_output2) <- paste("Iteration.",nResamples_seq,sep="")
coef_test_output2 <- data.frame(Observed.Values=val.output[,"LMA_gDW_m2"], 
                                coef_test_output2)

write.csv(coef_test_output2,file=file.path(outdir,paste0(inVar,'_Jackkife_PLSR_Val_Data_Output.csv')),
          row.names=FALSE)

## where does predint and confint come from?
write.csv(predint,file=file.path(outdir,paste0(inVar,'_PLSR_Val_Prediction_Intervals.csv')),
          row.names=FALSE)

write.csv(confint,file=file.path(outdir,paste0(inVar,'_PLSR_Val_Confidence_Intervals.csv')),
          row.names=FALSE)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### EOF