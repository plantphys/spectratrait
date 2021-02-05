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


#--------------------------------------------------------------------------------------------------#
### Load libraries
# make sure required tools are available 
req.packages <- c("devtools")
new.packages <- req.packages[!(req.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=c("Depends", "Imports",
                                                                       "LinkingTo"))
# install spectratrait package
devtools::install_github(repo = "TESTgroup-BNL/PLSR_for_plant_trait_prediction", dependencies=TRUE)
list.of.packages <- c("pls","dplyr","reshape2","here","plotrix","ggplot2","gridExtra","spectratrait")
invisible(lapply(list.of.packages, library, character.only = TRUE))
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Setup other functions and options
# not in
`%notin%` <- Negate(`%in%`)

# What is the source dataset from EcoSIS?
ecosis_id <- "960dbb0c-144e-4563-8117-9e23d14f4aa9"

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
### Example dataset
# 
# URL:
# https://ecosis.org/package/ngee-arctic-2016-leaf-spectral-reflectance-kougarok-road-watershed-seward-peninsula-alaska
#
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Get source dataset from EcoSIS
dat_raw <- spectratrait::get_ecosis_data(ecosis_id = ecosis_id)
head(dat_raw)
names(dat_raw)[1:40]
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
names(sample_info)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Plot spectra
spectratrait::f.plot.spec(Z=spectra,wv=wv,plot_label="NGEE-Arctic Leaf Spectra")
dev.copy(png,file.path(outdir,'Leaf_Spectra.png'), 
         height=2500,width=4900, res=340)
dev.off();
par(mfrow=c(1,1))
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Plot trait data
print("**** Plotting Ecosis trait data. Writing to scratch space ****")

# Organize leaf trait data
trait_data <- sample_info %>%
  select(Site,Sample_ID,USDA_Species_Code=`USDA Symbol`,
         Measurement_Date=`Sample Collection Date`,Cmass_g_g,Nmass_g_g,CN_Ratio,
         LMA_g_m2)
head(trait_data)        

# Prepare data for ggplot
trait_data <- melt(data = trait_data, id.vars = "USDA_Species_Code", measure.vars = c("LMA_g_m2",
                                                                                      "Cmass_g_g",
                                                                                      "Nmass_g_g",
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