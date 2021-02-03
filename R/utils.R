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
  handle <- textConnection(content(request, as = 'text'))
  on.exit(close(handle))
  read.table(handle, sep = sep, header = header)
}

##' Function to check for installed package
##' not presently used
testForPackage <- function(pkg) {
  if (!requireNamespace(pkg)) {
    stop("Package", pkg, "required but not installed")
  }
}