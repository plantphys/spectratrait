##' Create a calibration (training) / validation data split for PLSR model fitting and testing
##' 
##' @param dataset input full PLSR dataset to split into cal/val datasets
##' @param approach approach to splitting the dataset. Options: base or dplyr
##' @param split_seed random seed to use for splitting data
##' @param prop the proportion of data to preserve for calibration (e.g. 0.8) and validation (0.2). 
##' This sets the calibration proportion
##' @param group_variables Use factor variables to conduct a stratfied sampling for cal/val
##' 
##' @examples 
##' \dontrun{
##' print()
##' }
##' @return output_list A list containing the calibration dataset (cal_data)
##' and validation dataset (val_data)
##' 
##' @author Julien Lamour, Jeremiah Anderson, Shawn P. Serbin
##' @export
create_data_split <- function(dataset=NULL, approach=NULL, split_seed=123456789, prop=0.8,
                              group_variables=NULL) {
  set.seed(split_seed)
  if(!is.null(approach)) {
    if (approach=="base") {
      dataset$CalVal <- NA
      split_var <- group_variables
      if(length(group_variables) > 1){
        dataset$ID <- apply(dataset[, group_variables], MARGIN = 1, FUN = function(x) paste(x, collapse = " "))
      } else {
        dataset$ID <- dataset[, group_variables]
      }
      split_var_list <- unique(dataset$ID)
      for(i in 1:length(split_var_list)){
        temp <- row.names(dataset[ dataset$ID == split_var_list[i], ])
        ## there should probably be more than 4 obs I'm guessing, so this may need adjusting
        if(length(temp) > 3){
          Cal <- sample(temp,round(prop*length(temp)))
          Val <- temp[!temp %in% Cal]
          dataset$CalVal[ row.names(dataset) %in% Cal ] <- "Cal"
          dataset$CalVal[ row.names(dataset) %in% Val ] <- "Val"
          p_cal <- length(Cal)/length(temp) * 100
          message(paste0(split_var_list[i], "   ", "Cal", ": ", p_cal, "%"))
        } else {
          message(paste(split_var_list[i], "Not enough observations"))
        }
      }
      dataset$ID <- NULL
      # drop NA's in CalVal
      dataset <- dataset[!is.na(dataset$CalVal), ]
      cal.plsr.data <- dataset[dataset$CalVal== "Cal",]
      val.plsr.data <- dataset[dataset$CalVal== "Val",]
    } else if (approach=="dplyr")
      cal.plsr.data <- dataset %>% 
        group_by_at(vars(all_of(group_variables))) %>% 
        slice(sample(1:n(), prop*n())) %>% 
        data.frame()
    val.plsr.data <-dataset[!row.names(dataset) %in% row.names(cal.plsr.data),]
    
  } else {
    stop("**** Please choose either base R or dplyr data split ****")
  }
  output_list <- list(cal_data=cal.plsr.data, val_data=val.plsr.data)
  return(output_list)
}