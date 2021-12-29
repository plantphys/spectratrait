##' Function to source text data from GitHub
##' 
##' @param url http/https URL to the github dataset
##' @param sep dataset file delimiter
##' @param header TRUE/FALSE does the file have a column header?
##' 
##' @import httr
##' 
##' @author gist.github.com/christophergandrud/4466237
##' @export
source_GitHubData <- function(url, sep = ",", header = TRUE) {
  # define function to grab PLSR model from GitHub
  #devtools::source_gist("gist.github.com/christophergandrud/4466237")
  request <- httr::GET(url)
  httr::stop_for_status(request)
  handle <- textConnection(httr::content(request, as = 'text'))
  on.exit(close(handle))
  read.table(handle, sep = sep, header = header)
}

##' Calculate RMSE and percent RMSE with PLSR model results
##' 
##' @param plsr_dataset input plsr dataset
##' @param inVar the trait variable used in the calculation of RMSE
##' @param residuals predicted minus observed residual vector from either 
##' a cross-validation CV or independent validation
##' @param range calculate over the full data range or the 95% of data range. 
##' options full or 95perc
##' @return output a list containing the rmse and perc_rmse. 
##' output <- list(rmse = rmse, perc_rmse = perc_rmse)
##' 
##' @author Shawn P. Serbin
##' 
##' @export
percent_rmse <- function(plsr_dataset = NULL, inVar = NULL, 
                         residuals = NULL, range = "full") {
  rmse <- sqrt(mean(residuals^2, na.rm = TRUE))
  val_data_range <- range(plsr_dataset[,inVar], na.rm = TRUE)
  val_95perc_data_range <- quantile(plsr_dataset[,inVar], 
                                    probs = c(0.025,0.975), 
                                    na.rm = TRUE)
  if (range=="full") {
    perc_rmse <- (rmse/(val_data_range[2]-val_data_range[1]))*100
  } else if (range=="95perc") {
    perc_rmse <- (rmse/(val_95perc_data_range[[2]]-val_95perc_data_range[[1]]))*100 
  } else {
    print("Not a valid option, defaulting to full")
    perc_rmse <- (rmse/(val_data_range[2]-val_data_range[1]))*100
  }
  
  # output rmse list
  output <- list(rmse = rmse, perc_rmse = perc_rmse)
  return(output)
}


##' Not %in% function
##'
##' @export
`%notin%` <- Negate(`%in%`)

##' Function to check for installed package
##' not presently used
testForPackage <- function(pkg) {
  if (!requireNamespace(pkg)) {
    stop("Package", pkg, "required but not installed")
  }
}