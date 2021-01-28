####################################################################################################
#
#  
#   A simple example showing how to pull spectra and trait data from EcoSIS and plot it up
#
#   Spectra and trait data source:
#   https://ecosis.org/package/ngee-arctic-2016-leaf-spectral-reflectance-kougarok-road-seward-peninsula-alaska-2016
#
#    Notes:
#    * Provided as a basic example of how to apply the model to new spectra observations
#    * The author notes the code is not the most elegant or clean, but is functional 
#    * Questions, comments, or concerns can be sent to sserbin@bnl.gov
#    * Code is provided under GNU General Public License v3.0 
#
#
####################################################################################################


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files

list.of.packages <- c("readr","httr","dplyr","reshape2","ggplot2")  # packages needed for script
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=c("Depends", "Imports",
                                                                       "LinkingTo"))
# Load libraries
invisible(lapply(list.of.packages, library, character.only = TRUE))

#devtools::source_gist("gist.github.com/christophergandrud/4466237")
source_GitHubData <- function(url, sep = ",", header = TRUE) {
  require(httr)
  request <- GET(url)
  stop_for_status(request)
  handle <- textConnection(content(request, as = 'text'))
  on.exit(close(handle))
  read.table(handle, sep = sep, header = header)
}

# not in
`%notin%` <- Negate(`%in%`)

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
### Example datasets
# 
# URL:  https://ecosis.org/package/ngee-arctic-2016-leaf-spectral-reflectance-kougarok-road-seward-peninsula-alaska-2016
#
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Grab data
print("**** Downloading Ecosis data ****")
ecosis_id <- "960dbb0c-144e-4563-8117-9e23d14f4aa9"  # NGEE-Arctic dataset
ecosis_file <- sprintf(
  "https://ecosis.org/api/package/%s/export?metadata=true",
  ecosis_id
)
message("Downloading data...")
dat_raw <- read_csv(ecosis_file)
message("Download complete!")
names(dat_raw)[1:40]
head(dat_raw)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Prepare data
Start.wave <- 500
End.wave <- 2400
wv <- seq(Start.wave,End.wave,1)
spectra <- data.frame(dat_raw[,names(dat_raw) %in% wv])
names(spectra) <- c(paste0("Wave_",wv))
head(spectra)[,1:5]

sample_info <- dat_raw[,names(dat_raw) %notin% seq(350,2500,1)]
head(sample_info)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Plot spectra
cexaxis <- 1.5
cexlab <- 1.8
ylim <- 65
ylim2 <- 65

mean_spec <- colMeans(spectra[,which(names(spectra) %in% paste0("Wave_",wv))])
spectra_quantiles <- apply(spectra[,which(names(spectra) %in% paste0("Wave_",wv))],
                           2,quantile,na.rm=T,probs=c(0,0.025,0.05,0.5,0.95,0.975,1))
par(mfrow=c(1,1), mar=c(4.5,5.7,0.3,0.4), oma=c(0.3,0.9,0.3,0.1)) # B, L, T, R
plot(wv,mean_spec,ylim=c(0,ylim),cex=0.00001, col="white",xlab="Wavelength (nm)",
     ylab="Reflectance (%)",cex.axis=cexaxis, cex.lab=cexlab)
polygon(c(wv ,rev(wv)),c(spectra_quantiles[5,], rev(spectra_quantiles[3,])),
        col="#99CC99",border=NA)
lines(wv,mean_spec,lwd=3, lty=1, col="black")
lines(wv,spectra_quantiles[1,],lwd=1.85, lty=3, col="grey40")
lines(wv,spectra_quantiles[7,],lwd=1.85, lty=3, col="grey40")
legend("topright",legend=c("Mean reflectance","Min/Max", "95% CI"),lty=c(1,3,1),
       lwd=c(3,3,15),col=c("black","grey40","#99CC99"),bty="n", cex=1.7)
box(lwd=2.2)
dev.copy(png,file.path(outdir,
                       "NGEE-Arctic_2016_Kougarok_leaf_spectra_summary_plot.png"), 
         height=3000,width=3900, res=340)
dev.off();
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Plot trait data
print("**** Plotting Ecosis trait data. Writing to scratch space ****")

# Organize leaf trait data
names(sample_info)
trait_data <- sample_info %>%
  select(Site,Sample_ID,USDA_Species_Code=`USDA Symbol`,Common_Name=`Common Name`,LMA_gDW_m2=LMA_g_m2,
         Nmass_g_g,N_area_g_m2,Cmass_g_g,C_area_g_m2,CN_Ratio)
head(trait_data)

# Prepare data for ggplot
trait_data <- melt(data = trait_data, id.vars = "USDA_Species_Code", measure.vars = c("LMA_gDW_m2",
                                                                                      "Nmass_g_g",
                                                                                      "N_area_g_m2",
                                                                                      "CN_Ratio"))
head(trait_data)

# Graph the trait data and save a file to the scratch space
p2 <- ggplot(trait_data, aes(x=USDA_Species_Code, y=value)) + 
  geom_boxplot() +
  facet_wrap(~variable, scale="free")
print(p2)
ggsave(filename = file.path(outdir,"NGEE-Arctic_2016_Kougarok_Trait_data.png"), plot = p2,
       width = 40, height = 20, units = "cm")
#--------------------------------------------------------------------------------------------------#

### EOF