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
list.of.packages <- c("readr","httr","pls","dplyr","reshape2","ggplot2","gridExtra")  # packages needed for script
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load libraries
library(pls)
library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Setup other functions and options

# define function to grab PLSR model from GitHub
#devtools::source_gist("gist.github.com/christophergandrud/4466237")
source_GitHubData <-function(url, sep = ",", header = TRUE) {
  require(httr)
  request <- GET(url)
  stop_for_status(request)
  handle <- textConnection(content(request, as = 'text'))
  on.exit(close(handle))
  read.table(handle, sep = sep, header = header)
}

f.plot.spec=function(
  Z,                  ## Spectra matrix with each row corresponding to a spectra and wavelength in columns
  wv,                 ## vector of wavelengths corresponding to the column of the spectra matrix Z
  xlim=NULL,          ## vector to change the default xlim of the plots (ex xlim = c(500, 2400))
  position='topright',## Position of the legend (see base function legend for help)
  type='Reflectance', ## Name of the y axis and of the legend
  plot_label=NULL     ## optional label for plot
){
  if(mean(as.matrix(Z),na.rm=TRUE)>1){Z=Z/100} ## Check if the spectra are in pc [0,100] or in [0,1]
  if(is.null(xlim)){xlim=c(min(wv),max(wv))}
  mean_spec <- colMeans(Z)
  spectra_quantiles <- apply(Z,2,quantile,na.rm=T,probs=c(0,0.025,0.05,0.5,0.95,0.975,1))
  
  plot(x=NULL,y=NULL,ylim=c(0,100),xlim=xlim,xlab="Wavelength (nm)",
       ylab=type,main=plot_label)
  polygon(c(wv ,rev(wv)),c(spectra_quantiles[5,]*100, rev(spectra_quantiles[3,]*100)),
          col="#99CC99",border=NA)
  lines(wv,mean_spec*100,lwd=2, lty=1, col="black")
  lines(wv,spectra_quantiles[1,]*100, lty=3, col="grey40")
  lines(wv,spectra_quantiles[7,]*100, lty=3, col="grey40")
  legend(position,legend=c(paste("Mean",type),"Min/Max", "95% CI"),lty=c(1,3,1),
         lwd=c(2,1,10),col=c("black","grey40","#99CC99"),bty="n")
}

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
spectra <- data.frame(dat_raw[,names(dat_raw) %in% wv])
names(spectra) <- c(paste0("Wave_",wv))
head(spectra)[1:6,1:10]
sample_info <- dat_raw[,names(dat_raw) %notin% seq(350,2500,1)]
head(sample_info)

sample_info2 <- sample_info %>%
  select(Domain,Functional_type,Sample_ID,USDA_Species_Code=`USDA Symbol`,LMA_gDW_m2=LMA)
head(sample_info2)

plsr_data <- data.frame(sample_info2,spectra)
head(plsr_data)[,1:10]
rm(sample_info,sample_info2,spectra)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Create cal/val datasets
set.seed(2356812)

unique(plsr_data$USDA_Species_Code)
unique(plsr_data$Domain)
# !!! this is messy and could likely be streamlined !!!
# !!! also we may want to split data by both domain and functional type or species !!!
domains <- unique(plsr_data$Domain)
cal.plsr.data <- 0
val.plsr.data <- 0
prop <- 0.80
j <- 1
for (i in domains){
  print(paste("Domain: ",i,sep=""))
  temp.data <- plsr_data[which(plsr_data$Domain==i),]
  rows <- sample(1:nrow(temp.data),floor(prop*nrow(temp.data)))
  cal_data = droplevels(temp.data[rows,])
  val_data = droplevels(temp.data[-rows,])
  
  if(j==1){
    cal.plsr.data <- cal_data
    val.plsr.data <- val_data
  } else {
    cal.plsr.data <- rbind(cal.plsr.data,cal_data)
    val.plsr.data <- rbind(val.plsr.data,val_data)
  }
  
  j <- j+1
}
rm(temp.data)

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
full_plsr_data <- rbind(cal.plsr.data,val.plsr.data)
write.csv(full_plsr_data,file=paste0(outdir,inVar,'_Full_PLSR_Dataset.csv',sep=""),row.names=FALSE)
write.csv(cal.plsr.data,file=paste0(outdir,inVar,'_Cal_PLSR_Dataset.csv',sep=""),row.names=FALSE)
write.csv(val.plsr.data,file=paste0(outdir,inVar,'_Val_PLSR_Dataset.csv',sep=""),row.names=FALSE)
rm(cal_data,val_data,i,j,prop,rows,domains)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Create calibration and validation PLSR datasets

spec_start <- which(names(cal.plsr.data)==paste0("Wave_",Start.wave))
cal.spec <- as.matrix(droplevels(cal.plsr.data[,spec_start:dim(cal.plsr.data)[2]]))
cal.plsr.data.2 <- data.frame(cal.plsr.data[,1:spec_start-1],Spectra=I(cal.spec))
cal.plsr.data <- cal.plsr.data.2
head(cal.plsr.data)[,1:5]
rm(cal.plsr.data.2,cal.spec,spec_start)

spec_start <- which(names(val.plsr.data)==paste0("Wave_",Start.wave))
val.spec <- as.matrix(droplevels(val.plsr.data[,spec_start:dim(val.plsr.data)[2]]))
val.plsr.data.2 <- data.frame(val.plsr.data[,1:spec_start-1],Spectra=I(val.spec))
val.plsr.data <- val.plsr.data.2
head(val.plsr.data)[,1:5]
rm(val.plsr.data.2,val.spec,spec_start)

# plot cal and val spectra
par(mfrow=c(1,2)) # B, L, T, R
f.plot.spec(Z=cal.plsr.data$Spectra,wv=seq(Start.wave,End.wave,1),plot_label="Calibration")
f.plot.spec(Z=val.plsr.data$Spectra,wv=seq(Start.wave,End.wave,1),plot_label="Validation")
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Use Jackknife permutation to determine optimal number of compoenents
if(grepl("Windows", sessionInfo()$running)){
  pls.options(parallel =NULL)
} else {
  pls.options(parallel = parallel::detectCores()-1)
}
dims <- dim(plsr_data)
nComps <- 16
iterations <- 20
seg <- 15
prop <- 0.70
jk.out <- matrix(data=NA,nrow=iterations,ncol=nComps) 
print("*** Running jacknife permutation test.  Please hang tight, this can take awhile ***")
start.time <- Sys.time()
for (i in 1:iterations) {
  rows <- sample(1:nrow(cal.plsr.data),floor(prop*nrow(cal.plsr.data)))
  sub.data <- cal.plsr.data[rows,]
  plsr.out <- plsr(as.formula(paste(inVar,"~","Spectra")), scale=FALSE, center=TRUE, ncomp=nComps, 
                   validation="CV", segments = seg, segment.type="interleaved", trace=FALSE, data=sub.data)
  resPRESS <- as.vector(plsr.out$validation$PRESS)
  jk.out[i,seq(plsr.out$validation$ncomp)]=resPRESS
}
end.time <- Sys.time()
end.time - start.time

# Jackknife PRESS plot
pressDF <- as.data.frame(jk.out)
names(pressDF) <- as.character(seq(nComps))
pressDFres <- melt(pressDF)
bp <- ggplot(pressDFres, aes(x=variable, y=value)) + theme_bw() + 
  geom_boxplot(notch=TRUE) + labs(x="Number of Components", y="PRESS") +
  theme(axis.text=element_text(size=18), legend.position="none",
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))
bp

# conduct t.test across components to identify first minimum - just one of the ways to do this
j <-2 
results <- as.vector(array(data="NA", dim=c(nComps-1,1)))
for (i in seq_along(1:nComps-1)) {
  comp1 <- i; comp2 <- j
  ttest <- t.test(pressDFres$value[which(pressDFres$variable==comp1)],
                  pressDFres$value[which(pressDFres$variable==comp2)])
  #print(i)
  results[i] <- round(unlist(ttest$p.value),8)
  j <- j+1
  if (j > nComps) {
    break
  }
}
results <- data.frame(seq(2,nComps,1),results)
names(results) <- c("Component", "P.value")
results





#--------------------------------------------------------------------------------------------------#



