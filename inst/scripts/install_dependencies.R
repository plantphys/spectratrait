####################################################################################################
# Install dependencies
#
####################################################################################################


#--------------------------------------------------------------------------------------------------#
req.packages <- c("devtools","remotes","readr","RCurl","httr","pls","dplyr","reshape2","here",
                  "plotrix","scales","ggplot2","gridExtra")
new.packages <- req.packages[!(req.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=TRUE)
#--------------------------------------------------------------------------------------------------#