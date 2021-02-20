##' Run a PLSR model permutation analysis. Can be used to determine the optimal number of components 
##' or conduct a boostrap uncertainty analysis
##' 
##' See Serbin et al. (2019). DOI: https://doi.org/10.1111/nph.16123
##'
##' @param dataset input full PLSR dataset. Usually just the calibration dataset
##' @param maxComps maximum number of components to use for each PLSR fit
##' @param iterations how many different permutations to run
##' @param prop proportion of data to preserve for each permutation
##' @param verbose Should the function report the current iteration number to the terminal 
##' or run silently? TRUE/FALSE
##' 
##' @author Julien Lamour, Shawn P. Serbin
##' @export
pls_permutation <- function(dataset=NULL, maxComps=20, iterations=20, prop=0.70, 
                            verbose=FALSE) {
  coefs <- array(0,dim=c((ncol(dataset$Spectra)+1),iterations,maxComps))
  press.out <- array(data=NA, dim=c(iterations,maxComps))
  print("*** Running permutation test.  Please hang tight, this can take awhile ***")
  print("Options:")
  print(paste("Max Components:",maxComps, "Iterations:", iterations, 
              "Data Proportion:", prop, sep=" "))
  # TODO: ADD PROGRESS BAR HERE!!!
  for (i in seq_along(1:iterations)) {
    if (verbose) {
      message(paste("Running interation", i))
    }
    rows <- sample(1:nrow(dataset),floor(prop*nrow(dataset)))
    sub.data <- dataset[rows,]
    val.sub.data <- dataset[-rows,]
    plsr.out <- plsr(as.formula(paste(inVar,"~","Spectra")), scale=FALSE, center=TRUE, 
                     ncomp=maxComps, validation="none", data=sub.data)
    pred_val <- predict(plsr.out,newdata=val.sub.data)
    sq_resid <- (pred_val[,,]-val.sub.data[,inVar])^2
    press <- apply(X = sq_resid, MARGIN = 2, FUN = sum)
    press.out[i,] <- press
    coefs[,i,] <- coef(plsr.out, intercept = TRUE, ncomp = 1:maxComps)
    rm(rows,sub.data,val.sub.data,plsr.out,pred_val,sq_resid,press)
  }
  # create a new list with PRESS and permuted coefficients x wavelength x component number
  output <- list(PRESS=press.out, coef_array=coefs)
  return(output)
}