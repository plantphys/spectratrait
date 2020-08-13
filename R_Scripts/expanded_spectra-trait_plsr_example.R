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
lapply(list.of.packages, require,character.only = TRUE)
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
  "https://ecosis.org/api/package/%s/export?metadata=true",
  ecosis_id
)
message("Downloading data...")
dat_raw <- read_csv(ecosis_file)
message("Download complete!")
head(dat_raw)
names(dat_raw)[1:40]
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Create full plsr dataset
Start.wave <- 500
End.wave <- 2400
wv <- seq(Start.wave,End.wave,1)
## Majuscul to be consitent everywhere
Spectra <- as.matrix(dat_raw[,names(dat_raw) %in% wv])
colnames(Spectra) <- c(paste0("Wave_",wv))
head(Spectra)[1:6,1:10]
sample_info <- dat_raw[,names(dat_raw) %notin% seq(350,2500,1)] ## Or sample_info <- dat_raw[,!names(dat_raw) %in% seq(350,2500,1)] so you dont have to define %notin% earlier
head(sample_info)

sample_info2 <- sample_info %>%
  select(Domain,Functional_type,Sample_ID,USDA_Species_Code=`USDA Symbol`,LMA_gDW_m2=LMA)
head(sample_info2)

plsr_data <- data.frame(sample_info2,Spectra=I(Spectra))
rm(sample_info,sample_info2,Spectra)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Create cal/val datasets
set.seed(2356812)

## See the proportion of samples per species and Domain
table(plsr_data$USDA_Species_Code,plsr_data$Domain)

## Make a stratified random sampling in the strata USDA_Species_Code and Domain
prop=0.8
cal.plsr.data <- plsr_data %>% group_by(USDA_Species_Code,Domain) %>% slice(sample(1:n(), prop*n())) %>% data.frame()
val.plsr.data <- plsr_data[!plsr_data$Sample_ID %in% cal.plsr.data$Sample_ID,]

## Verification of the stratified sampling, are the proportion similar with the plsr_data dataset?
table(cal.plsr.data$USDA_Species_Code,cal.plsr.data$Domain)

# Datasets:
print(paste("Cal observations: ",dim(cal.plsr.data)[1],sep=""))
print(paste("Val observations: ",dim(val.plsr.data)[1],sep=""))

cal_hist_plot <- qplot(cal.plsr.data[,paste0(inVar)],geom="histogram",binwidth = 10,
      main = paste0("Cal. Histogram for ",inVar),
      xlab = paste0(inVar),ylab = "Count",fill=I("grey50"),col=I("black"),alpha=I(.7))
val_hist_plot <- qplot(val.plsr.data[,paste0(inVar)],geom="histogram",binwidth = 10,
      main = paste0("Val. Histogram for ",inVar),
      xlab = paste0(inVar),ylab = "Count",fill=I("grey50"),col=I("black"),alpha=I(.7))
grid.arrange(cal_hist_plot, val_hist_plot, ncol=2)

# !!  do we need to actually write any of this out to temp dir? !!
## For me the answer is no!
## I dont really like when a script write files on my computer. I feel that it is intrusive ^^

#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Create calibration and validation PLSR datasets
# Already done earlier?

# plot cal and val spectra
par(mfrow=c(1,2)) # B, L, T, R
f.plot.spec(Z=cal.plsr.data$Spectra,wv=seq(Start.wave,End.wave,1),plot_label="Calibration")
f.plot.spec(Z=val.plsr.data$Spectra,wv=seq(Start.wave,End.wave,1),plot_label="Validation")
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#

### Use permutation to determine the optimal number of components

if(grepl("Windows", sessionInfo()$running)){
  pls.options(parallel =NULL)
} else {
  pls.options(parallel = parallel::detectCores()-1)
}
dims <- dim(plsr_data)
nComps <- 20
iterations <- 20
seg <- 15

plsr.out=plsr(as.formula(paste(inVar,"~","Spectra")), scale=FALSE, center=TRUE, ncomp=nComps, 
     validation="CV", segments = seg, segment.type="interleaved", trace=FALSE, data=cal.plsr.data)

nComps=selectNcomp(plsr.out, method = "onesigma", plot = TRUE)

fit <- plsr.out$fitted.values[,1,nComps]
pls.options(parallel = NULL)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
#ALready done for using selectNcomp


#calibration
cal_plot_data <- data.frame(cal.plsr.data[, which(names(cal.plsr.data) %notin% "Spectra")], Fitted=fit)
cal_plot_data <- cal_plot_data %>%
  mutate(Residuals = Fitted-LMA_gDW_m2)
head(cal_plot_data)

# validation
val_plot_data <- data.frame(val.plsr.data[, which(names(val.plsr.data) %notin% "Spectra")], 
                            Fitted=as.vector(predict(plsr.out, newdata = val.plsr.data, ncomp=nComps, 
                                                     type="response")[,,1]))
val_plot_data <- val_plot_data %>%
  mutate(Residuals = Fitted-LMA_gDW_m2)
head(val_plot_data)

# cal
cal_scatter_plot <- ggplot(cal_plot_data, aes(x=Fitted, y=LMA_gDW_m2)) + 
  theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1, color="dark grey", 
                                          linetype="dashed", size=1.5) + xlim(0, 275) + ylim(0, 275) +
  labs(x=expression(paste("Predicted LMA (",g~m^{-2},")")), 
       y=expression(paste("Observed LMA (",g~m^{-2},")"))) +
  annotate("text", x=250, y=70, label = paste0("R^2 == ", round(pls::R2(plsr.out)[[1]][nComps],2)), parse=T) + 
  annotate("text", x=250, y=40, label = paste0("RMSE == ", round(pls::RMSEP(plsr.out)[[1]][nComps],2)), parse=T) +
  annotate("text",x=20,y=220,label="Calibration") + 
  theme(axis.text=element_text(size=18), legend.position="none",
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))

cal_resid_histogram <- ggplot(cal_plot_data, aes(x=Residuals)) +
  geom_histogram(binwidth=.5, alpha=.5, position="identity") + 
  geom_vline(xintercept = 0, color="black", 
             linetype="dashed", size=1) + theme_bw() + 
  theme(axis.text=element_text(size=18), legend.position="none",
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))

# val
val_scatter_plot <- ggplot(val_plot_data, aes(x=Fitted, y=LMA_gDW_m2)) + 
  theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1, color="dark grey", 
                                          linetype="dashed", size=1.5) + xlim(0, 275) + ylim(0, 275) +
  labs(x=expression(paste("Predicted LMA (",g~m^{-2},")")), 
       y=expression(paste("Observed LMA (",g~m^{-2},")"))) +
  annotate("text", x=250, y=70, label = paste0("R^2 == ", 
                                               round(pls::R2(plsr.out, newdata = val.plsr.data)[[1]][nComps],2)), parse=T) + 
  annotate("text", x=250, y=40, label = paste0("RMSE == ", 
                                               round(pls::RMSEP(plsr.out, newdata = val.plsr.data)[[1]][nComps],2)), parse=T) +
  annotate("text",x=20,y=220,label="Validation") + 
  theme(axis.text=element_text(size=18), legend.position="none",
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))

val_resid_histogram <- ggplot(val_plot_data, aes(x=Residuals)) +
  geom_histogram(binwidth=.5, alpha=.5, position="identity") + 
  geom_vline(xintercept = 0, color="black", 
             linetype="dashed", size=1) + theme_bw() + 
  theme(axis.text=element_text(size=18), legend.position="none",
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))

# plot cal/val side-by-side
grid.arrange(cal_scatter_plot, val_scatter_plot, cal_resid_histogram, val_resid_histogram, 
             nrow=2,ncol=2)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### results by functional type and domain
# validation by functional type
scatter_plot <- ggplot(val_plot_data, aes(x=Fitted, y=LMA_gDW_m2)) + 
  theme_bw() + geom_point(aes(fill=Functional_type),alpha=0.6,colour="black", pch=21, size=4) + 
  geom_abline(intercept = 0, slope = 1, color="dark grey", 
              linetype="dashed", size=1.5) + xlim(0, 275) + ylim(0, 275) +
  labs(x=expression(paste("Predicted LMA (",g~m^{-2},")")), 
       y=expression(paste("Observed LMA (",g~m^{-2},")"))) +
  annotate("text", x=250, y=70, label = paste0("R^2 == ", 
                                               round(pls::R2(plsr.out, newdata = val.plsr.data)[[1]][nComps],2)), 
           parse=T) + 
  annotate("text", x=250, y=40, label = paste0("RMSE == ", 
                                               round(pls::RMSEP(plsr.out, newdata = val.plsr.data)[[1]][nComps],2)), 
           parse=T) +
  theme(axis.text=element_text(size=18), legend.position="bottom",legend.title=element_text(size=16),
        legend.text=element_text(size=14),
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))
scatter_plot

resid_histogram <- ggplot(val_plot_data, aes(x=Residuals, fill=Functional_type)) +
  geom_histogram(binwidth=.5, alpha=.5, position="identity") + 
  geom_vline(xintercept = 0, color="black", alpha=0.6,
             linetype="dashed", size=1) + theme_bw() + 
  theme(axis.text=element_text(size=18), legend.position="bottom",legend.title=element_text(size=16),
        legend.text=element_text(size=14),
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))
resid_histogram

# NEON domain
scatter_plot <- ggplot(val_plot_data, aes(x=Fitted, y=LMA_gDW_m2)) + 
  theme_bw() + geom_point(aes(fill=Domain),alpha=0.6,colour="black", pch=21, size=4) + 
  geom_abline(intercept = 0, slope = 1, color="dark grey", 
              linetype="dashed", size=1.5) + xlim(0, 275) + ylim(0, 275) +
  labs(x=expression(paste("Predicted LMA (",g~m^{-2},")")), 
       y=expression(paste("Observed LMA (",g~m^{-2},")"))) +
  annotate("text", x=250, y=70, label = paste0("R^2 == ", 
                                               round(pls::R2(plsr.out, newdata = val.plsr.data)[[1]][nComps],2)), 
           parse=T) + 
  annotate("text", x=250, y=40, label = paste0("RMSE == ", 
                                               round(pls::RMSEP(plsr.out, newdata = val.plsr.data)[[1]][nComps],2)), 
           parse=T) +
  theme(axis.text=element_text(size=18), legend.position="bottom",legend.title=element_text(size=16),
        legend.text=element_text(size=14),
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))
scatter_plot

resid_histogram <- ggplot(val_plot_data, aes(x=Residuals, fill=Domain)) +
  geom_histogram(binwidth=.5, alpha=.5, position="identity") + 
  geom_vline(xintercept = 0, color="black", alpha=0.6,
             linetype="dashed", size=1) + theme_bw() + 
  theme(axis.text=element_text(size=18), legend.position="bottom",legend.title=element_text(size=16),
        legend.text=element_text(size=14),
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))
resid_histogram
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
dev.off()
par(mfrow=c(2,1))
plot(seq(Start.wave,End.wave,1),coefs,cex=0.01,xlab="Wavelength (nm)",ylab="REG COEF")
lines(seq(Start.wave,End.wave,1),coefs,lwd=2.5)
abline(h=0,lty=2,col="dark grey")

plot(seq(Start.wave,End.wave,1),vips,xlab="Wavelength (nm)",ylab="VIP",cex=0.01)
lines(seq(Start.wave,End.wave,1),vips,lwd=3)
abline(h=0.8,lty=2,col="dark grey")
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
coefs <- coef(plsr.out,ncomp=nComps,intercept=TRUE)
write.csv(coefs,file=file.path(outdir,paste0(inVar,'_PLSR_Coefficients_',nComps,'comp.csv')),
          row.names=TRUE)

# PLSR VIP
write.csv(vips,file=file.path(outdir,paste0(inVar,'_PLSR_VIPs_',nComps,'comp.csv')))

# confirm files were written to temp space
print("**** PLSR output files: ")
list.files(getwd())[grep(pattern = inVar, list.files(getwd()))]
#--------------------------------------------------------------------------------------------------#


#---------------- Jackknife model evaluation ------------------------------------------------------#
#!!  this code section needs lots of cleaning and refining.  not optimal !!
nComps
resamples <- 100 #1000
output.jackknife.stats <- data.frame(Rsq=rep(NA,resamples),RMSEP=rep(NA,resamples),
                                     PERC_RMSEP=rep(NA,resamples), Bias=rep(NA,resamples))
output.jackknife.coefs <- array(data=NA,dim=c(resamples,
                                              dim(coef(plsr.out,ncomp=nComps,intercept=TRUE))[1]))
output.jackknife.coefs.scaled <- array(data=NA,dim=c(resamples,
                                                     dim(coef(plsr.out,ncomp=nComps,intercept=TRUE))[1]))
vips <- array(data=NA,dim=c(resamples,
                            dim(coef(plsr.out,ncomp=nComps,intercept=FALSE))[1]))

for (i in 1:resamples) {
  rows <- sample(1:nrow(cal.plsr.data),floor(0.7*nrow(cal.plsr.data)))
  cal.data.jk <- cal.plsr.data[rows,]
  val.data.jk <- cal.plsr.data[-rows,]
  
  dimsCal <- dim(cal.data.jk)
  dimsVal <- dim(val.data.jk)
  
  ### Build PLSR model with training data
  if(grepl("Windows", sessionInfo()$running)){
    pls.options(parallel =NULL)
  } else {
    pls.options(parallel = parallel::detectCores()-1)
  }
  pls.jack <- plsr(as.formula(paste0(inVar,"~","Spectra")), scale=FALSE, ncomp=nComps, validation="none", data=cal.data.jk)
  pls.jack.scaled <- plsr(as.formula(paste0(inVar,"~","Spectra")), scale=TRUE, ncomp=nComps, validation="none", data=cal.data.jk)
  
  ### Estimate independent (Validation) samples
  plsr.val <- val.data.jk[,inVar]
  pred.val.data <- as.vector(predict(pls.jack,newdata=val.data.jk$Spectra,
                                      ncomp=nComps,type="response")[,,1])
  
  ### Coefficients and VIPs
  output.jackknife.coefs[i,] <- as.vector(coef(pls.jack,ncomp=nComps,intercept=TRUE))
  output.jackknife.coefs.scaled[i,] <- as.vector(coef(pls.jack.scaled,ncomp=nComps,intercept=TRUE))
  vips[i,] <- VIP(pls.jack)[nComps,]
  
  # Error statistics
  n <- length(plsr.val)
  Rsq.val <- R2(pls.jack,newdata = val.data.jk)$val[,,nComps+1] 
  val.residuals <- pred.val.data-plsr.val
  Val.bias <- mean(pred.val.data)-mean(plsr.val)
  MSEP <- mean(val.residuals^2)
  RMSEP.val <- sqrt(MSEP)
  PERC_RMSEP <- (RMSEP.val/(max(plsr.val)-min(plsr.val)))*100
  
  ### Store results of iteration i
  output.jackknife.stats[i,1] <- Rsq.val
  output.jackknife.stats[i,2] <- RMSEP.val
  output.jackknife.stats[i,3] <- PERC_RMSEP
  output.jackknife.stats[i,4] <- Val.bias
  
  print(paste("Running Iteration",i))
  print(paste("Stats: ","Rsq ",round(Rsq.val,2)," / RMSEP ",round(RMSEP.val,2), " / %RMSEP ",
              round(PERC_RMSEP,2)," / Bias ",round(Val.bias,2),sep="" ) )
  flush.console()  # force the output
  
  # Remove temp objects
  rm(cal.data.jk,val.data.jk,n,pls.jack,plsr.val,Rsq.val,pred.val.data,RMSEP.val,PERC_RMSEP,
     val.residuals,Val.bias)
}
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Histogram statistics
head(output.jackknife.stats)
rsqHist <- qplot(output.jackknife.stats[,"Rsq"],geom="histogram",
                 main = "Jackknife Rsq",
                 xlab = paste0(inVar),ylab = "Count",fill=I("grey50"),col=I("black"),alpha=I(.7))
rmseHist <- qplot(output.jackknife.stats[,"RMSEP"],geom="histogram",
                 main = "Jackknife RMSE",
                 xlab = paste0(inVar),ylab = "Count",fill=I("grey50"),col=I("black"),alpha=I(.7))
biasHist <- qplot(output.jackknife.stats[,"Bias"],geom="histogram",
                  main = "Jackknife Bias",
                  xlab = paste0(inVar),ylab = "Count",fill=I("grey50"),col=I("black"),alpha=I(.7))
grid.arrange(rsqHist, rmseHist, biasHist, ncol=3)
#--------------------------------------------------------------------------------------------------#


#---------------- Output jackknife results --------------------------------------------------------#
write.csv(output.jackknife.stats,file=file.path(outdir,paste0(inVar,'_Jackkife_PLSR_Resutls.csv')),
          row.names=FALSE)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### plot jackknife coefficients
dev.off()
dims <- dim(output.jackknife.coefs)
plot.coefs <- output.jackknife.coefs[,2:dims[2]]
plot.min <- min(output.jackknife.coefs[,2:dims[2]])
plot.max <- max(output.jackknife.coefs[,2:dims[2]])
jk.intercepts <- output.jackknife.coefs[,1]

# Stats
coef.means <- colMeans(plot.coefs)
sd.coef <- apply(plot.coefs,MARGIN=2,FUN=function(x)sd(x))
min.coef <- apply(plot.coefs,MARGIN=2,FUN=function(x)min(x))
max.coef <- apply(plot.coefs,MARGIN=2,FUN=function(x)max(x))
coef.quant <- apply(plot.coefs,2,quantile,probs=c(0.025,0.975))
intercepts.quant <- quantile(jk.intercepts,probs=c(0.025,0.975))

# T Test
x <- resamples # From Jackknife above
results <- apply(plot.coefs,2, function(plot.coefs) {
  t.test(x = plot.coefs[1:x])$p.value}) 

waves <- seq(Start.wave,End.wave,1)
coefs <- as.vector(coef(plsr.out,ncomp=nComps,intercept=FALSE))
plot(waves,coefs,type="l",lwd=4,ylim=c(plot.min,plot.max))
# Min/Max
polygon(c(waves ,rev(waves)),c(max.coef, rev(min.coef)),col="grey50",border=NA)
lines(waves,min.coef,lty=1,lwd=3,col="grey50")
lines(waves,max.coef,lty=1,lwd=3,col="grey50")
# 95% CIs
polygon(c(waves ,rev(waves)),c(coef.quant[2,], rev(coef.quant[1,])),col="grey70",border=NA)
lines(waves,coef.quant[1,],lty=1,lwd=2,col="grey70")
lines(waves,coef.quant[2,],lty=1,lwd=2,col="grey70")
# replot the mean and zero line
lines(waves,coefs,lwd=4)
abline(h=0,lty=2,col="grey",lwd=1.5)
legend("bottomright",legend=c("Mean","95% CI"),lty=c(1,1),
       col=c("black","dark grey"),lwd=3)
box(lwd=2.2)

# P-vals
plot(waves,results,pch=21,bg="grey80",ylab="P-value",xlab="Wavelength (nm)",
     cex=2)
#--------------------------------------------------------------------------------------------------#


#---------------- Output jackknife results --------------------------------------------------------#
# JK Coefficents
out.jk.coefs <- data.frame(Iteration=seq(1,resamples,1),jk.intercepts,plot.coefs)
names(out.jk.coefs) <- c("Iteration","Intercept",paste("Wave_",seq(Start.wave,End.wave,1),sep=""))
write.csv(out.jk.coefs,file=file.path(outdir,paste0(inVar,'_Jackkife_PLSR_Coefficients.csv')),
          row.names=FALSE)

# VIPs
out.jk.vips <- data.frame(Iteration=seq(1,resamples,1),vips)
names(out.jk.vips) <- c("Iteration",paste("Wave_",seq(Start.wave,End.wave,1),sep=""))
write.csv(out.jk.vips,file=file.path(outdir,paste0(inVar,'_Jackkife_PLSR_VIPs.csv')),
          row.names=FALSE)

# Coeff quantiles
out.coef.quant <- array(data=NA,dim=c(2,dim(out.jk.coefs)[2]))
out.coef.quant[1,1] <- "5%"
out.coef.quant[2,1] <- "95%"
out.coef.quant[1,2] <- intercepts.quant[[1]]
out.coef.quant[2,2] <- intercepts.quant[[2]]
out.coef.quant[,3:dim(out.jk.coefs)[2]] <- coef.quant
out.coef.quant <- data.frame(out.coef.quant)

names(out.coef.quant) <- c("Quantile","Intercept",paste("Wave_",seq(Start.wave,End.wave,1),sep=""))

write.csv(out.coef.quant,file=file.path(outdir,paste0(inVar,'_Jackkife_PLSR_Coefficient_Quantiles.csv')),
          row.names=TRUE)
# P-vals
out.pvals <- data.frame(Wavelength=paste("Wave_",seq(Start.wave,End.wave,1),sep=""),Pval=results)
write.csv(out.pvals,file=file.path(outdir,paste0(inVar,'_Jackkife_PLSR_Coefficient_Pvals.csv')),
          row.names=FALSE)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# JK Val plot
dims <- dim(output.jackknife.coefs)
intercepts <- output.jackknife.coefs[,1]
jk.coef.test.output <- array(data=NA,dim=c(dim(val.plsr.data)[1],dims[1]))
for (i in 1:length(intercepts)){
  coefs <- as.vector(output.jackknife.coefs[i,2:dims[2]])
  temp <- val.plsr.data$Spectra %*% coefs  # Updated: Using matrix mult.
  vals <- data.frame(rowSums(temp))+intercepts[i]
  jk.coef.test.output[,i] <- vals[,1]
}
pred.quant <- apply(jk.coef.test.output,1,quantile,probs=c(0.025,0.975))
pred.quant.ll <- pred.quant[1,]
pred.quant.ul <- pred.quant[2,]

jk_val_plot_data <- data.frame(val.output, LL=pred.quant.ll, UL=pred.quant.ul)
head(jk_val_plot_data)

# plot needs cleaning up. draft
jk_val_scatterplot <- ggplot(jk_val_plot_data, aes(x=PLSR_Predicted, y=LMA_gDW_m2)) + 
  theme_bw() + geom_point() + geom_errorbar(aes(xmin = LL,xmax = UL), width = 0.2) + 
  geom_abline(intercept = 0, slope = 1, color="dark grey", 
            linetype="dashed", size=1.5) + xlim(0, 275) + ylim(0, 275) +
  labs(x=expression(paste("Predicted LMA (",g~m^{-2},")")), 
       y=expression(paste("Observed LMA (",g~m^{-2},")"))) + 
  theme(axis.text=element_text(size=18), legend.position="none",
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))
jk_val_scatterplot
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Output JK Coefficient test results
jk.coef.test.output2 <- data.frame(jk.coef.test.output)
names(jk.coef.test.output2) <- paste("Iteration.",seq(1,resamples,1),sep="")
jk.coef.test.output2 <- data.frame(Observed.Values=val.output[,"LMA_gDW_m2"], 
                                   jk.coef.test.output2)

write.csv(jk.coef.test.output2,file=file.path(outdir,paste0(inVar,'_Jackkife_PLSR_Val_Data_Output.csv')),
          row.names=FALSE)

write.csv(predint,file=file.path(outdir,paste0(inVar,'_PLSR_Val_Prediction_Intervals.csv')),
          row.names=FALSE)

write.csv(confint,file=file.path(outdir,paste0(inVar,'_PLSR_Val_Confidence_Intervals.csv')),
          row.names=FALSE)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### EOF