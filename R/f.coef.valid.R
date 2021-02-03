## Return the intercept and the coefficients of the jackknife validation
##' @title f.coef.valid
##' 
##' @param plsr.out plsr model obtained with jaccknife = TRUE 
##' @param data_plsr data used for the plsr model with Spectra the matrix of spectra
##' @param ncomp number of selection components
##' @param inVar Name of the PLSR model response variable
##' 
##' @return B returns the intercept and the coefficients of the jackknife or bootstrap validation 
##' 
##' @author Julien Lamour
##' @export
f.coef.valid <- function(plsr.out, data_plsr, ncomp, inVar) {
  ## Only work in the case where center=TRUE in the plsr model
  B <- plsr.out$validation$coefficients[, , ncomp,, drop = FALSE]
  dB <- dim(B)
  dB[1] <- dB[1] + 1
  dnB <- dimnames(B)
  dnB[[1]] <- c("(Intercept)", dnB[[1]])
  BInt <- array(dim = dB, dimnames = dnB)
  BInt[-1, , ,] <- B
  nseg=dB[[4]]
  for (i in 1:nseg){
    Y<-data_plsr[,inVar]
    Y<-Y[-plsr.out$validation$segments[[i]]]
    Ymeans<-mean(Y)
    X<-data_plsr$Spectra
    X<-X[-plsr.out$validation$segments[[i]],]
    Xmeans<-colMeans(X)
    BInt[1, , ,i] <- Ymeans - Xmeans %*% B[, , , i]
  } 
  B <- BInt
  return(B)
}