####################################################################################################
#
#  
#   A simple example "How-to" script illustrating the use of basic PLSR modeling to develop a 
#   spectra-trait algorithm to estimate leaf mass area with leaf-level spectroscopy data. The 
#   example is built from published data source from the EcoSIS spectral database. 
#
#   Spectra and trait data source:
#   https://ecosis.org/package/fresh-leaf-spectra-to-estimate-lma-over-neon-domains-in-eastern-united-states
#
#    Notes:
#    * Provided as a basic example of how to apply the model to new spectra observations
#    * The author notes the code is not the most elegant or clean, but is functional 
#    * Questions, comments, or concerns can be sent to sserbin@bnl.gov
#    * Code is provided under GNU General Public License v3.0 
#
#
#    --- Last updated:  11.19.2019 By Shawn P. Serbin <sserbin@bnl.gov>
####################################################################################################


## !!! THIS IS A VERY BASIC, PRELIM SCRIPT.  NEEDS WORK TO GENERALIZE ON TESTING ON MULTIPLE PLATFORMS

#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files

list.of.packages <- c("devtools","readr","scales","plotrix","httr","pls","dplyr",
                      "reshape2")  # packages needed for script
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load libraries needed for script
library(readr)    # readr - read_csv function to pull data from EcoSIS
library(plotrix)  # plotCI - to generate obsvered vs predicted plot with CIs
library(scales)   # alpha() - for applying a transparency to data points
library(devtools)
library(pls)
library(dplyr)
library(reshape2)

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

# not in
`%notin%` <- Negate(`%in%`)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Set working directory (scratch space)
output_dir <- file.path("~",'scratch/')
if (! file.exists(output_dir)) dir.create(output_dir,recursive=TRUE, showWarnings = FALSE)
setwd(output_dir) # set working directory
getwd()  # check wd
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Example datasets
# 
# URL:  https://ecosis.org/package/fresh-leaf-spectra-to-estimate-lma-over-neon-domains-in-eastern-united-states
#
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Grab data
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
## Create PLSR dataset
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

plsr_data <- data.frame(sample_info2,Spectra=I(as.matrix(spectra)))
str(plsr_data)

inVar <- "LMA_gDW_m2"
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Jackknife test to find number of components - this could be optimized
dims <- dim(plsr_data)
nComps <- 15
iterations <- 50
seg <- 15
prop <- 0.70
jk.out <- matrix(data=NA,nrow=iterations,ncol=nComps) 
start.time <- Sys.time()
for (i in 1:iterations) {
  print(paste("Iteration: ",i,sep=""))
  rows <- sample(1:nrow(plsr_data),floor(prop*nrow(plsr_data)))
  sub.data <- plsr_data[rows,]
  plsr.out <- plsr(as.formula(paste(inVar,"~","Spectra")), scale=FALSE, ncomp=nComps, validation="CV",
                   segments = seg, segment.type="random", trace=15, data=sub.data)
  resPRESS <- as.vector(plsr.out$validation$PRESS)
  jk.out[i,seq(plsr.out$validation$ncomp)]=resPRESS
}
end.time <- Sys.time()
end.time - start.time

pressDF <- as.data.frame(jk.out)
names(pressDF) <- as.character(seq(nComps))
pressDFres <- melt(pressDF)
png(paste(output_dir,inVar,"_findComponents.png",sep=""),width=6,height=5,units="in",res=200)
boxplot(pressDFres$value~pressDFres$variable,xlab="n Components",ylab="PRESS",main=inVar)
dev.off()

# How many components? - !!automate this!!
loc.1 <- 14
loc.2 <- 15
ttest <- t.test(pressDFres$value[which(pressDFres$variable==loc.1)],pressDFres$value[which(pressDFres$variable==loc.2)])
ttest

# Test components
nComps <- 14
#plsr.out <- plsr(as.formula(paste(inVar,"~","Spectra")),scale=FALSE,ncomp=nComps,validation="LOO",
#                 trace=TRUE,data=cal.plsr.data)
segs <- 100
plsr.out <- plsr(as.formula(paste(inVar,"~","Spectra")),scale=FALSE,ncomp=nComps,validation="CV",
                 segments=segs, segment.type="random",trace=15,data=plsr_data)
fit1 <- plsr.out$fitted.values[,1,nComps]
#temp <- cal.plsr.data[,inVar]
#temp <- temp[complete.cases(temp)]
plot(as.vector(fit1),plsr_data[,inVar],xlim=c(0,range(plsr_data[,inVar])[2]+1),
     ylim=c(0,range(plsr_data[,inVar])[2]+1))
abline(lm(plsr_data[,inVar]~fit1),lwd=2)
abline(0,1,col="red",lwd=2,lty=2)
#--------------------------------------------------------------------------------------------------#






