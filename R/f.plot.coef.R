##' @title f.plot.coef
##' 
##' @param Z Coefficient matrix with each row corresponding to the coefficients and wavelength in columns
##' @param wv vector of wavelengths
##' @param xlim vector to change the default xlim of the plots (ex xlim = c(500, 2400))
##' @param position Position of the legend (see base function legend for help)
##' @param type Name of the y axis and of the legend
##' @param plot_label optional plot label to include with the figure
##' 
##' @importFrom stats quantile
##' @importFrom graphics polygon lines legend box
##' 
##' @author Julien Lamour
##' @export
f.plot.coef <- function(
  Z,                  ## Coefficient matrix with each row corresponding to the coefficients and wavelength in columns
  wv,                 ## vector of wavelengths 
  xlim=NULL,          ## vector to change the default xlim of the plots (ex xlim = c(500, 2400))
  position='topright',## Position of the legend (see base function legend for help)
  type='Coefficient', ## Name of the y axis and of the legend
  plot_label=NULL     ## optional label for plot
){
  
  if(is.null(xlim)){xlim=c(min(wv),max(wv))}
  mean_spec <- colMeans(Z)
  spectra_quantiles <- apply(Z,2,quantile,na.rm=T,probs=c(0,0.01,0.025,0.05,0.5,0.95,0.975,0.99,1))
  
  plot(x=NULL,y=NULL,xlim=xlim,ylim=c(min(Z),max(Z)),xlab="Wavelength (nm)",
       ylab=type,main=plot_label)
  polygon(c(wv ,rev(wv)),c(spectra_quantiles[9,], rev(spectra_quantiles[1,])),
          col="grey60",border=NA)
  polygon(c(wv ,rev(wv)),c(spectra_quantiles[6,], rev(spectra_quantiles[4,])),
          col="#99CC99",border=NA)
  lines(wv,mean_spec,lwd=2, lty=1, col="black")
  lines(wv,spectra_quantiles[1,], lty=3, col="grey60")
  lines(wv,spectra_quantiles[9,], lty=3, col="grey60")
  legend(position,legend=c(paste("Mean",type),"Min/Max (range)", "95% CI"),lty=c(1,1,1),
         lwd=c(2,10,10),col=c("black","grey50","#99CC99"),bty="n")
  box(lwd=2.2)
}