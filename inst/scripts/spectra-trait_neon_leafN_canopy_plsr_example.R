# canopy example

#--------------------------------------------------------------------------------------------------#
### Step 1.
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
### Step 2.
### Setup other functions and options
# not in
`%notin%` <- Negate(`%in%`)

# Script options
pls::pls.options(plsralg = "oscorespls")
pls::pls.options("plsralg")

# Default par options
opar <- par(no.readonly = T)

# What is the target variable?
inVar <- "Nitrogen"

# What is the source dataset from EcoSIS?
ecosis_id <- "b9dbf3db-5b9c-4ab2-88c2-26c8b39d0903"

# Specify output directory, output_dir 
# Options: 
# tempdir - use a OS-specified temporary directory 
# user defined PATH - e.g. "~/scratch/PLSR"
output_dir <- "tempdir"
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Step 3.
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
### Step 4.
### Get source dataset from EcoSIS
dat_raw <- spectratrait::get_ecosis_data(ecosis_id = ecosis_id)
head(dat_raw)
names(dat_raw)[1:40]
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Step 5.
### Create plsr dataset

# identify the trait data and other metadata
sample_info <- dat_raw[,names(dat_raw) %notin% seq(300,2600,1)]
head(sample_info)

# spectra matrix
Spectra <- as.matrix(dat_raw[,names(dat_raw) %notin% names(sample_info)])
head(Spectra)[1:6,1:10]

# set the desired spectra wavelength range to include
Start.wave <- 500
End.wave <- 2400
wv <- seq(Start.wave,End.wave,1)
final_spec <- Spectra[,round(as.numeric(colnames(Spectra))) %in% wv]
colnames(final_spec) <- c(paste0("Wave_",colnames(final_spec)))
head(final_spec)[1:6,1:10]

# Drop bad data
dropwaves <- c(1360:1440, 1826:1946)
final_spec <- final_spec[,colnames(final_spec) %notin% paste0("Wave_",dropwaves)]
head(final_spec)[1:6,1:10]
wv <- as.numeric(gsub(pattern = "Wave_",replacement = "", x = colnames(final_spec)))

# assemble dataset !! ADDD MORE HERE !!
sample_info2 <- sample_info %>%
  select(Plot_ID,Sample_Year,SLA,Nitrogen)

##!! need to split PLOT_ID to retain first part for selecting samples

plsr_data <- data.frame(sample_info2,final_spec)
rm(sample_info,sample_info2,Spectra)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Step 6.
#### Example data cleaning.  End user needs to do what's appropriate for their 
#### data.  This may be an iterative process.
# Keep only complete rows of inVar and spec data before fitting
plsr_data <- plsr_data[complete.cases(plsr_data[,names(plsr_data) %in% 
                                                  c(inVar,paste0("Wave_",wv))]),]
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Step 7.
### Create cal/val datasets
## Make a stratified random sampling in the strata USDA_Species_Code and Domain

method <- "dplyr" #base/dplyr
# base R - a bit slow
# dplyr - much faster
split_data <- spectratrait::create_data_split(dataset=plsr_data, approach=method, split_seed=7529075,
                                              prop=0.8, group_variables="Species_Code")
names(split_data)
cal.plsr.data <- split_data$cal_data
head(cal.plsr.data)[1:8]
val.plsr.data <- split_data$val_data
head(val.plsr.data)[1:8]
rm(split_data)


