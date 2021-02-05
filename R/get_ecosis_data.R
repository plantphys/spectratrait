##' Function to pull data from EcoSIS using the EcoSIS API
##' 
##' @param ecosis_id the alphanumeric EcoSIS API dataset ID
##' 
##' @examples 
##' \dontrun{ 
##' ecosis_id <- "960dbb0c-144e-4563-8117-9e23d14f4aa9"
##' dat_raw <- get_ecosis_data(ecosis_id = ecosis_id)
##' head(dat_raw)
##' names(dat_raw)[1:40]
##' }
##' 
##' @return EcoSIS spectral dataset object
##' 
##' @author Shawn P. Serbin, Alexey Shiklomanov
##' @export
get_ecosis_data <- function(ecosis_id = NULL) {
  if(!is.null(ecosis_id)) {
    print("**** Downloading Ecosis data ****")
    ecosis_id <- ecosis_id
    ecosis_file <- sprintf(
      "https://ecosis.org/api/package/%s/export?metadata=true",
      ecosis_id)
    message("Downloading data...")
    dat_raw <- readr::read_csv(ecosis_file)
    message("Download complete!")
    return(dat_raw)
  } else {
    stop("**** No EcoSIS ID provided.  Please provide a valid ID before proceeding ****")
  }
}