##Same function as pls_permutation but just with an additional argument (groups),
## and a new sampling approach which takes into account the input groups

##Stratifying is only included if groups is not null
##A dplyr approach is used to do the by-group sampling


## groups should be an vector of caracters of the form c("var1", "var2"..."varn")
##I have not included any control for this condition, though

##Changed commands marked 



pls_permut_by_groups <- function (dataset = NULL, maxComps = 20, iterations = 20, prop = 0.7, 
          verbose = FALSE, groups=NULL) 
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
##New condition statement to select the sampling strategy (random or stratified) based on
##the content of the group argument
    if (!is.null(groups)) {
##A new dataset is created for internal training by
##creating a new Internal Id at the beginning,
##grouping by the groups object,
##and sampling within each group a number of cases equal to prop*n()
##Then the selected values of the internal id is stored as the rows vector,
##which is used to separate the training dataset and the validation data set
##(unchanged from pls_permutation)
      trainset<- dataset %>%
        mutate(int_id = row_number()) %>%
        group_by_at(groups) %>%
        slice(sample(1:n(), prop * n())) 
      rows<- trainset$int_id
    } else {
        rows <- sample(1:nrow(dataset), floor(prop * nrow(dataset)))
    }
    sub.data <- dataset[rows, ]
    val.sub.data <- dataset[-rows, ]
    plsr.out <- plsr(as.formula(paste(inVar, "~", "Spectra")), 
                     scale = FALSE, center = TRUE, ncomp = maxComps, validation = "none", 
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
