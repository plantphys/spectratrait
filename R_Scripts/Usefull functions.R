

f.plot.spec=function(
  Z,                  ## Spectra matrix with each row corresponding to a spectra and wavelenght in columns
  wv,                 ## vector of wavelengths corresponding to the column of the spectra matrix Z
  xlim=NULL,          ## vector to change the default xlim of the plots (ex xlim = c(500, 2400))
  position='topright',## Position of the legend (see base function legend for help)
  type='Reflectance'  ## Name of the y axis and of the legend
){
  if(mean(as.matrix(Z),na.rm=TRUE)>1){Z=Z/100} ## Check if the spectra are in pc [0,100] or in [0,1]
  if(is.null(xlim)){xlim=c(min(wv),max(wv))}
  mean_spec <- colMeans(Z)
  spectra_quantiles <- apply(Z,2,quantile,na.rm=T,probs=c(0,0.025,0.05,0.5,0.95,0.975,1))
  
  plot(x=NULL,y=NULL,ylim=c(0,100),xlim=xlim,xlab="Wavelength (nm)",
       ylab=type)
  polygon(c(wv ,rev(wv)),c(spectra_quantiles[5,]*100, rev(spectra_quantiles[3,]*100)),
          col="#99CC99",border=NA)
  lines(wv,mean_spec*100,lwd=2, lty=1, col="black")
  lines(wv,spectra_quantiles[1,]*100, lty=3, col="grey40")
  lines(wv,spectra_quantiles[7,]*100, lty=3, col="grey40")
  legend(position,legend=c(paste("Mean",type),"Min/Max", "95% CI"),lty=c(1,3,1),
         lwd=c(2,1,10),col=c("black","grey40","#99CC99"),bty="n")
}
