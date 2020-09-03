####################################################################################################
#
#
#  Helper functions
#
#  Last updated: 2020-07-29 by Shawn P. Serbin
####################################################################################################


#--------------------------------------------------------------------------------------------------#
# define function to grab PLSR model from GitHub
#devtools::source_gist("gist.github.com/christophergandrud/4466237")
source_GitHubData <-function(url, sep = ",", header = TRUE) {
  require(httr)
  request <- GET(url)
  stop_for_status(request)
  handle <- textConnection(content(request, as = 'text'))
  on.exit(close(handle))
  read.table(handle, sep = sep, header = header)
}
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
#' @title f.plot.spec
#' @param Z Spectra matrix with each row corresponding to a spectra and wavelength in columns
#' @param wv vector of wavelengths corresponding to the column of the spectra matrix Z
#' @param xlim vector to change the default xlim of the plots (ex xlim = c(500, 2400))
#' @param position Position of the legend (see base function legend for help)
#' @param type Name of the y axis and of the legend. E.g. Reflectance, Transmittance
#' @param plot_label option plot label
#' @author Julien Lamour
#' 
f.plot.spec=function(
  Z,                  ## Spectra matrix with each row corresponding to a spectra and wavelength in columns
  wv,                 ## vector of wavelengths corresponding to the column of the spectra matrix Z
  xlim=NULL,          ## vector to change the default xlim of the plots (ex xlim = c(500, 2400))
  position='topright',## Position of the legend (see base function legend for help)
  type='Reflectance', ## Name of the y axis and of the legend
  plot_label=NULL     ## optional label for plot
){
  if(mean(as.matrix(Z),na.rm=TRUE)>1){Z=Z/100} ## Check if the spectra are in pc [0,100] or in [0,1]
  if(is.null(xlim)){xlim=c(min(wv),max(wv))}
  mean_spec <- colMeans(Z)
  spectra_quantiles <- apply(Z,2,quantile,na.rm=T,probs=c(0,0.025,0.05,0.5,0.95,0.975,1))
  
  plot(x=NULL,y=NULL,ylim=c(0,100),xlim=xlim,xlab="Wavelength (nm)",
       ylab=type,main=plot_label)
  polygon(c(wv ,rev(wv)),c(spectra_quantiles[5,]*100, rev(spectra_quantiles[3,]*100)),
          col="#99CC99",border=NA)
  lines(wv,mean_spec*100,lwd=2, lty=1, col="black")
  lines(wv,spectra_quantiles[1,]*100, lty=3, col="grey40")
  lines(wv,spectra_quantiles[7,]*100, lty=3, col="grey40")
  legend(position,legend=c(paste("Mean",type),"Min/Max", "95% CI"),lty=c(1,3,1),
         lwd=c(2,1,10),col=c("black","grey40","#99CC99"),bty="n")
}
#--------------------------------------------------------------------------------------------------#
f.plot.coef=function(
  Z,                  ## Spectra matrix with each row corresponding to a spectra and wavelength in columns
  wv,                 ## vector of wavelengths corresponding to the column of the spectra matrix Z
  xlim=NULL,          ## vector to change the default xlim of the plots (ex xlim = c(500, 2400))
  position='topright',## Position of the legend (see base function legend for help)
  type='Reflectance', ## Name of the y axis and of the legend
  plot_label=NULL     ## optional label for plot
){
  
  if(is.null(xlim)){xlim=c(min(wv),max(wv))}
  mean_spec <- colMeans(Z)
  spectra_quantiles <- apply(Z,2,quantile,na.rm=T,probs=c(0,0.025,0.05,0.5,0.95,0.975,1))
  
  plot(x=NULL,y=NULL,xlim=xlim,ylim=c(min(Z),max(Z)),xlab="Wavelength (nm)",
       ylab=type,main=plot_label)
  polygon(c(wv ,rev(wv)),c(spectra_quantiles[5,], rev(spectra_quantiles[3,])),
          col="#99CC99",border=NA)
  lines(wv,mean_spec,lwd=2, lty=1, col="black")
  lines(wv,spectra_quantiles[1,], lty=3, col="grey40")
  lines(wv,spectra_quantiles[7,], lty=3, col="grey40")
  legend(position,legend=c(paste("Mean",type),"Min/Max", "95% CI"),lty=c(1,3,1),
         lwd=c(2,1,10),col=c("black","grey40","#99CC99"),bty="n")
}

#--------------------------------------------------------------------------------------------------#
## VIP returns all VIP values for all variables and all number of components,
## as a ncomp x nvars matrix.
VIP <- function(object) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*") # Replace with matrix mult.
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## VIPjh returns the VIP of variable j with h components
VIPjh <- function(object, j, h) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  b <- c(object$Yloadings)[1:h]
  T <- object$scores[,1:h, drop = FALSE]
  SS <- b^2 * colSums(T^2)
  W <- object$loading.weights[,1:h, drop = FALSE]
  Wnorm2 <- colSums(W^2)
  sqrt(nrow(W) * sum(SS * W[j,]^2 / Wnorm2) / sum(SS))
}
#--------------------------------------------------------------------------------------------------#