## code to prepare `DATASET` dataset goes here
#usethis::use_data(DATASET, overwrite = TRUE)
getwd()
devtools::load_all()
library(dplyr)
`%notin%` <- Negate(`%in%`)

# get data
ecosis_id <- "25770ad9-d47c-428b-bf99-d1543a4b0ec9"
dat_raw <- spectratrait::get_ecosis_data(ecosis_id = ecosis_id)
head(dat_raw)
names(dat_raw)[1:60]

# setup data
Start.wave <- 500
End.wave <- 2400
wv <- seq(Start.wave,End.wave,1)
Spectra <- as.matrix(dat_raw[,names(dat_raw) %in% wv])
colnames(Spectra) <- c(paste0("Wave_",wv))
head(Spectra)[1:6,1:10]
sample_info <- dat_raw[,names(dat_raw) %notin% seq(350,2500,1)]
head(sample_info)

sample_info2 <- sample_info %>%
  select(Species_Code=`USDA Symbol`, Common_Name=`Common Name`, C_N_mass, C_g_m2, 
         H20_g_m2, LMA_g_m2, N_g_m2, QC)

plsr_data <- data.frame(sample_info2,Spectra)
plsr_data <- plsr_data %>%
  filter(is.na(QC)) %>%
  select(-QC)
plsr_data <- plsr_data[complete.cases(plsr_data),]

write.csv(x = plsr_data, file = "ely_plsr_data.csv", row.names = F)
ely_plsr_data <- read.csv("ely_plsr_data.csv", header = T)

usethis::use_data(ely_plsr_data, overwrite = TRUE, internal = FALSE)
                              
                              
                              