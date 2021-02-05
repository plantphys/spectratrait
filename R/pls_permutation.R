##' Run a PLSR model permutation analysis. Can be used to determine the optimal number of components 
##' or conduct a boostrap uncertainty analysis
##' 
##' See Serbin et al. (2019). DOI: https://doi.org/10.1111/nph.16123
##'
##' @param dataset input full PLSR dataset. Usually just the calibration dataset
##' @param maxComps maximum number of components to use for each PLSR fit
##' @param iterations how many different permutations to run
##' @param seg currently unused - should be removed from this function call
##' @param prop proportion of data to preserve for each permutation
##' 
##' @author Julien Lamour, Shawn P. Serbin
##' @export
pls_permutation <- function(dataset=NULL, maxComps=20, iterations=20, seg=100, prop=0.70) {
  coefs <- array(0,dim=c((ncol(dataset$Spectra)+1),iterations,maxComps))
  press.out <- array(data=NA, dim=c(iterations,maxComps))
  print("*** Running permutation test.  Please hang tight, this can take awhile ***")
  print(paste("Options:", maxComps, iterations, seg, prop, sep=" "))
  for (i in seq_along(1:iterations)) {
    message(paste("Running interation", i))
    rows <- sample(1:nrow(dataset),floor(prop*nrow(dataset)))
    sub.data <- dataset[rows,]
    val.sub.data <- dataset[-rows,]
    plsr.out <- plsr(as.formula(paste(inVar,"~","Spectra")), scale=FALSE, center=TRUE, ncomp=maxComps, 
                     validation="none", data=sub.data)
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