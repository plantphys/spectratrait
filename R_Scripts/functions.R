####################################################################################################
#
#
#  Helper functions
#
#
####################################################################################################


#--------------------------------------------------------------------------------------------------#
# define function to grab PLSR model from GitHub
#devtools::source_gist("gist.github.com/christophergandrud/4466237")
source_GitHubData <- function(url, sep = ",", header = TRUE) {
  require(httr)
  request <- GET(url)
  stop_for_status(request)
  handle <- textConnection(content(request, as = 'text'))
  on.exit(close(handle))
  read.table(handle, sep = sep, header = header)
}
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
#' @author Shawn P. Serbin
get_ecosis_data <- function(ecosis_id = NULL) {
  if(!is.null(ecosis_id)) {
    print("**** Downloading Ecosis data ****")
    ecosis_id <- ecosis_id
    ecosis_file <- sprintf(
      "https://ecosis.org/api/package/%s/export?metadata=true",
      ecosis_id)
    message("Downloading data...")
    dat_raw <- read_csv(ecosis_file)
    message("Download complete!")
    return(dat_raw)
  } else {
    stop("**** No EcoSIS ID provided.  Please provide a valid ID before proceeding ****")
  }
}
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
#' @author Julien Lamour, Jeremiah Anderson, Shawn P. Serbin
create_data_split <- function(approach=NULL, split_seed=123456789, prop=0.8,
                              group_variables=NULL) {
  set.seed(split_seed)
  if(!is.null(approach)) {
    if (approach=="base") {
      plsr_data$CalVal <- NA
      split_var <- group_variables
      if(length(group_variables) > 1){
        plsr_data$ID <- apply(plsr_data[, group_variables], MARGIN = 1, FUN = function(x) paste(x, collapse = " "))
      } else {
        plsr_data$ID <- plsr_data[, group_variables]
      }
      split_var_list <- unique(plsr_data$ID)
      for(i in 1:length(split_var_list)){
        temp <- row.names(plsr_data[ plsr_data$ID == split_var_list[i], ])
        ## there should probably be more than 4 obs I'm guessing, so this may need adjusting
        if(length(temp) > 3){
          Cal <- sample(temp,round(prop*length(temp)))
          Val <- temp[!temp %in% Cal]
          plsr_data$CalVal[ row.names(plsr_data) %in% Cal ] <- "Cal"
          plsr_data$CalVal[ row.names(plsr_data) %in% Val ] <- "Val"
          p_cal <- length(Cal)/length(temp) * 100
          message(paste0(split_var_list[i], "   ", "Cal", ": ", p_cal, "%"))
        } else {
          message(paste(split_var_list[i], "Not enough observations"))
        }
      }
      plsr_data$ID <- NULL
      # drop NA's in CalVal
      plsr_data <- plsr_data[!is.na(plsr_data$CalVal), ]
      cal.plsr.data <- plsr_data[plsr_data$CalVal== "Cal",]
      val.plsr.data <- plsr_data[plsr_data$CalVal== "Val",]
    } else if (approach=="dplyr")
      cal.plsr.data <- plsr_data %>% 
        group_by_at(vars(all_of(group_variables))) %>% 
        slice(sample(1:n(), prop*n())) %>% 
        data.frame()
      #val.plsr.data <- plsr_data[!plsr_data$Sample_ID %in% cal.plsr.data$Sample_ID,]
      val.plsr.data <-plsr_data[!row.names(plsr_data) %in% row.names(cal.plsr.data),]
    
  } else {
    stop("**** Please choose either base R or dplyr data split ****")
  }
  output_list <- list(cal_data=cal.plsr.data, val_data=val.plsr.data)
  return(output_list)
}
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
#' @title pls_permutation
#' @author Julien Lamour, Shawn P. Serbin
#'
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
    coefs[,i,] <- coef(plsr.out,intercept = TRUE, ncomp = 1:maxComps)
    rm(rows,sub.data,val.sub.data,plsr.out,pred_val,sq_resid,press)
  }
  # create a new list with PRESS and permuted coefficients x wavelength x component number
  output <- list(PRESS=press.out, coef_array=coefs)
  return(output)
}
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
#' @param method Which approach to use to find optimal components. Options: pls, firstPlateau, firstMin
#' 
#' @author Julien Lamour, Jeremiah Anderson, Shawn P. Serbin
#' 
find_optimal_components <- function(dataset=NULL, method="pls", maxComps=20, iterations=20, seg=100, 
                                    prop=0.70, random_seed=123456789) {
  set.seed(random_seed)
  if(method=="pls") {
    print("*** Running PLS permutation test ***")
    maxComps <- maxComps
    seg <- seg
    plsr.out <- pls::plsr(as.formula(paste(inVar,"~","Spectra")), scale=FALSE, center=TRUE, ncomp=maxComps, 
                     validation="CV", segments = seg, segment.type="interleaved", trace=FALSE, 
                     jackknife=TRUE, data=cal.plsr.data)
    nComps <- selectNcomp(plsr.out, method = "onesigma", plot = TRUE)
  }
  if(method=="firstPlateau") {
    press.out <- pls_permutation(dataset=dataset, maxComps=maxComps, iterations=iterations, 
                                 seg=seg, prop=prop)
    # PRESS plot
    pressDF <- as.data.frame(press.out$PRESS)
    names(pressDF) <- as.character(seq(maxComps))
    pressDFres <- reshape2::melt(pressDF)
    results <- NULL
    for(i in 1:(maxComps-1)){
      p_value <- t.test(press.out$PRESS[,i], press.out$PRESS[,(i+1)])$p.value
      temp_results <- data.frame(Component=(i+1), P.value= round(p_value, 6))
      results <- rbind(results, temp_results)
    }
    nComps <- min(results[results$P.value > 0.05, "Component"])
    print(paste0("*** Optimal number of components based on t.test: ", nComps))
    bp <- ggplot(pressDFres, aes(x=variable, y=value)) + theme_bw() + 
      geom_boxplot(notch=FALSE) + labs(x="Number of Components", y="PRESS") + 
      stat_boxplot(geom = "errorbar", width = 0.2) +
      geom_vline(xintercept = nComps, linetype="dashed", 
                 color = "blue", size=1.0)
    theme(axis.text=element_text(size=18), legend.position="none",
          axis.title=element_text(size=20, face="bold"), 
          axis.text.x = element_text(angle = 0,vjust = 0.5),
          panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))
    print(bp)
  }
  if(method=="firstMin") {
    press.out <- pls_permutation(dataset=dataset, maxComps=maxComps, iterations=iterations, 
                                 seg=seg, prop=prop)
    # PRESS plot
    pressDF <- as.data.frame(press.out$PRESS)
    names(pressDF) <- as.character(seq(maxComps))
    pressDFres <- reshape2::melt(pressDF)
    # find lowest press
    mean_PRESS_comp <- apply(X = pressDF, MARGIN = 2, FUN = mean)
    lowest_PRESS <- which.min(mean_PRESS_comp)
    results <- as.vector(array(data="NA", dim=c(lowest_PRESS-1,1)))
    for (i in seq_along(1:(lowest_PRESS-1))) {
      comp1 <- i; comp2 <- lowest_PRESS
      ttest <- t.test(pressDFres$value[which(pressDFres$variable==comp1)],
                      pressDFres$value[which(pressDFres$variable==comp2)])
      #print(i)
      results[i] <- round(unlist(ttest$p.value),8)
    }
    results <- data.frame(seq(1,lowest_PRESS-1,1),results)
    names(results) <- c("Component", "P.value")
    first <- min(which(as.numeric(as.character(results$P.value)) > 0.05))
    nComps <- results$Component[first]
    print(paste0("*** Optimal number of components based on t.test: ", nComps))
    bp <- ggplot(pressDFres, aes(x=variable, y=value)) + theme_bw() + 
      geom_boxplot(notch=FALSE) + labs(x="Number of Components", y="PRESS") + 
      stat_boxplot(geom = "errorbar", width = 0.2) +
      geom_vline(xintercept = nComps, linetype="dashed", 
                 color = "blue", size=1.0)
    theme(axis.text=element_text(size=18), legend.position="none",
          axis.title=element_text(size=20, face="bold"), 
          axis.text.x = element_text(angle = 0,vjust = 0.5),
          panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))
    print(bp)
  }
  return(nComps)
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
f.plot.spec <- function(
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
       ylab=paste0(type," (%)"),main=plot_label)
  polygon(c(wv ,rev(wv)),c(spectra_quantiles[5,]*100, rev(spectra_quantiles[3,]*100)),
          col="#99CC99",border=NA)
  lines(wv,mean_spec*100,lwd=2, lty=1, col="black")
  lines(wv,spectra_quantiles[1,]*100, lty=3, col="grey40")
  lines(wv,spectra_quantiles[7,]*100, lty=3, col="grey40")
  legend(position,legend=c(paste("Mean",type),"Min/Max", "95% CI"),lty=c(1,3,1),
         lwd=c(2,1,10),col=c("black","grey40","#99CC99"),bty="n")
}
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
#' @title f.plot.coef
#' @author Julien Lamour
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
  legend(position,legend=c(paste("Mean",type),"Min/Max", "95% CI"),lty=c(1,3,1),
         lwd=c(2,1,10),col=c("black","grey50","#99CC99"),bty="n")
}
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
## Return the intercept and the coefficients of the jackknife validation
## Only work in the case where center=TRUE in the plsr model
#' @author Julien Lamour
f.coef.valid <- function(plsr.out, ## plsr model with center = TRUE
                       data_plsr, ## data used for the plsr model
                       ncomp ## number of selection components
){
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
#--------------------------------------------------------------------------------------------------#

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

#--------------------------------------------------------------------------------------------------#
### EOF