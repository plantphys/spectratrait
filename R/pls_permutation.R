##' Run a PLSR model permutation analysis. Can be used to determine the optimal number of components 
##' or conduct a boostrap uncertainty analysis
##' 
##' See Serbin et al. (2019). DOI: https://doi.org/10.1111/nph.16123
##'
##' @param dataset input full PLSR dataset. Usually just the calibration dataset
##' @param maxComps maximum number of components to use for each PLSR fit
##' @param iterations how many different permutations to run
##' @param prop proportion of data to preserve for each permutation
##' @param verbose Should the function report the current iteration status/progress to the terminal
##' or run silently? TRUE/FALSE. Default FALSE
##' @return output a list containing the PRESS and coef_array.
##' output <- list(PRESS=press.out, coef_array=coefs)
##' 
##' @importFrom pls plsr 
##' @importFrom utils flush.console read.table setTxtProgressBar txtProgressBar
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
              "Data Proportion (percent):", prop*100, sep=" "))
  
  if (verbose) {
    j <- 1 # <--- Numeric counter for progress bar
    pb <- txtProgressBar(min = 0, max = iterations, 
                         char="*",width=70,style = 3)
  }
  
  for (i in seq_along(1:iterations)) {
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
    
    ### Display progress to console
    if (verbose) {
      setTxtProgressBar(pb, j)    # show progress bar
      j <- j+1                    # <--- increase counter by 1
      flush.console()             #<--- show output in real-time
    }
  }
  if (verbose) {
    close(pb)
  }
  
  # create a new list with PRESS and permuted coefficients x wavelength x component number
  print("*** Providing PRESS and coefficient array output ***")
  output <- list(PRESS=press.out, coef_array=coefs)
  return(output)
}



##Same function as pls_permutation but just with an additional argument (groups),
## and a new sampling approach which takes into account the input groups

##Stratifying is only included if groups is not null
##A dplyr approach is used to do the by-group sampling


## groups should be an vector of caracters of the form c("var1", "var2"..."varn")
##I have not included any control for this condition, though

##Changed commands marked 

##' Run a PLSR model permutation analysis stratified by selected "groups". Can be used to 
##' determine the optimal number of components or conduct a boostrap uncertainty analysis
##' 
##' @param dataset input full PLSR dataset. Usually just the calibration dataset
##' @param maxComps maximum number of components to use for each PLSR fit
##' @param iterations how many different permutations to run
##' @param prop proportion of data to preserve for each permutation
##' @param verbose Should the function report the current iteration status/progress to the terminal
##' or run silently? TRUE/FALSE. Default FALSE
##' @param group_variables Use factor variables to conduct a stratified sampling for the
##' PLSR permutation analysis. group_variables should be of the form c("var1", "var2"..."varn")
##' 
##' @return output a list containing the PRESS and coef_array.
##' output <- list(PRESS=press.out, coef_array=coefs)
##' 
##' @import dplyr
##' @importFrom pls plsr 
##' @importFrom utils flush.console read.table setTxtProgressBar txtProgressBar
##' 
##' @author asierrl, Shawn P. Serbin, Julien Lamour
##' @export
##' 
pls_permutation_by_groups <- function (dataset = NULL, maxComps = 20, iterations = 20, prop = 0.7, 
                                  verbose = FALSE, group_variables=NULL) 
{
  coefs <- array(0, dim = c((ncol(dataset$Spectra) + 1), iterations, 
                            maxComps))
  press.out <- array(data = NA, dim = c(iterations, maxComps))
  print("*** Running permutation test.  Please hang tight, this can take awhile ***")
  print("Options:")
  print(paste("Max Components:", maxComps, "Iterations:", 
              iterations, "Data Proportion (percent):", prop * 
                100, sep = " "))
  if (verbose) {
    j <- 1
    pb <- utils::txtProgressBar(min = 0, max = iterations, 
                                char = "*", width = 70, style = 3)
  }
  for (i in seq_along(1:iterations)) {
    if (!is.null(group_variables)) {
      trainset <- dataset %>%
        mutate(int_id = row_number()) %>%
        group_by_at(group_variables) %>%
        slice(sample(1:n(), prop * n()))
      rows <- trainset$int_id
      } else {
       rows <- sample(1:nrow(dataset), floor(prop * nrow(dataset)))
    }
    sub.data <- dataset[rows, ]
    val.sub.data <- dataset[-rows, ]
    plsr.out <- plsr(as.formula(paste(inVar, "~", "Spectra")), 
                     scale = FALSE, center = TRUE, ncomp = maxComps, 
                     validation = "none", 
                     data = sub.data)
    pred_val <- predict(plsr.out, newdata = val.sub.data)
    sq_resid <- (pred_val[, , ] - val.sub.data[, inVar])^2
    press <- apply(X = sq_resid, MARGIN = 2, FUN = sum)
    press.out[i, ] <- press
    coefs[, i, ] <- coef(plsr.out, intercept = TRUE, ncomp = 1:maxComps)
    rm(rows, sub.data, val.sub.data, plsr.out, pred_val, 
       sq_resid, press)
    if (verbose) {
      setTxtProgressBar(pb, j)
      j <- j + 1
      flush.console()
      }
  }
  if (verbose) {
    close(pb)
  }
  print("*** Providing PRESS and coefficient array output ***")
  output <- list(PRESS = press.out, coef_array = coefs)
  return(output)
}