##' Applies different methods for the determination of the optimal number of PLSR model components
##' 
##' @param dataset input full PLSR dataset. Usually just the calibration dataset
##' @param targetVariable What object or variable to use as the Y (predictand) in the PLSR model? 
##' Usually the "inVar" variable set at the beginning of a PLS script
##' @param method Which approach to use to find optimal components. Options: pls, firstPlateau, firstMin
##' @param maxComps maximum number of components to consider
##' @param iterations how many different permutations to run
##' @param seg For the built-in pls method, how many different data segments to select from the input dataset
##' @param prop proportion of data to preserve for each permutation
##' @param random_seed random seed to use for splitting data
##' 
##' @importFrom stats as.formula coef predict quantile t.test
##' @importFrom pls plsr selectNcomp
##' @importFrom reshape2 melt
##' @import ggplot2
##' 
##' @return nComps the optimal number of PLSR components
##' 
##' @author Julien Lamour, Jeremiah Anderson, Shawn P. Serbin
##' @export 
find_optimal_components <- function(dataset=NULL, targetVariable=NULL, method="pls", maxComps=20, 
                                    iterations=20, seg=100, prop=0.70, random_seed=123456789) {
  
  set.seed(random_seed)
  inVar <- targetVariable
  print("*** Identifying optimal number of PLSR components ***")
  
  if(method=="pls") {
    print("*** Running PLS permutation test ***")
    
    plsr.out <- pls::plsr(as.formula(paste(inVar,"~","Spectra")), scale=FALSE, center=TRUE, ncomp=maxComps, 
                          validation="CV", segments = seg, segment.type="interleaved", trace=FALSE, 
                          jackknife=TRUE, data=dataset)
    nComps <- pls::selectNcomp(plsr.out, method = "onesigma", plot = TRUE)
  }
  if(method=="firstPlateau") {
    press.out <- spectratrait::pls_permutation(dataset=dataset, targetVariable=inVar, maxComps=maxComps, 
                                               iterations=iterations, prop=prop)
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
    press.out <- spectratrait::pls_permutation(dataset=dataset, targetVariable=inVar, maxComps=maxComps, 
                                               iterations=iterations, prop=prop)
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

##' Uses the firstMin and firstPlateau methods for the determination of the optimal number of PLSR model components,
##' by group (i.e. optimal selection by stratification)
##' 
##' @param dataset input full PLSR dataset. Usually just the calibration dataset
##' @param targetVariable What object or variable to use as the Y (predictand) in the PLSR model? 
##' Usually the "inVar" variable set at the beginning of a PLS script
##' @param method Which approach to use to find optimal components. Options: firstPlateau, firstMin
##' @param maxComps maximum number of components to consider
##' @param iterations how many different permutations to run
##' @param prop proportion of data to preserve for each permutation
##' @param random_seed random seed to use for splitting data
##' @param group_variables group_variables character vector of the form c("var1", "var2"..."varn") 
##' providing the factors used for stratified sampling. 
##' 
##' @importFrom stats as.formula coef predict quantile t.test
##' @import ggplot2
##' @importFrom reshape2 melt
##' 
##' @return nComps the optimal number of PLSR components
##' 
##' @author asierrl, Shawn P. Serbin
##' @export 
find_optimal_comp_by_groups <- function (dataset = NULL, targetVariable = NULL, method = "firstPlateau",
                                         maxComps = 20, iterations = 20, prop = 0.7, random_seed = 123456789, 
                                         group_variables = NULL) {
  set.seed(random_seed)
  inVar <- targetVariable
  # TODO - really should merge this with the original and have an if/else if not NULL and select either
  # pls_permutation OR pls_permutation_by_groups.
  print("*** Identifying optimal number of PLSR components using stratified resampling by group_variables ***")
  if (method == "pls") {
    stop("*** Please select either the firstMin and firstPlateau. The pls package approach is not compatible ***")
  }
  if (method == "firstPlateau") {
    press.out <- spectratrait::pls_permutation_by_groups(dataset=dataset, targetVariable=inVar,
                                                    maxComps=maxComps, iterations=iterations, 
                                                    prop=prop, group_variables=group_variables)
    pressDF <- as.data.frame(press.out$PRESS)
    names(pressDF) <- as.character(seq(maxComps))
    pressDFres <- reshape2::melt(pressDF)
    results <- NULL
    for (i in 1:(maxComps - 1)) {
      p_value <- t.test(press.out$PRESS[, i], press.out$PRESS[, (i + 1)])$p.value
      temp_results <- data.frame(Component = (i + 1), P.value = round(p_value, 6))
      results <- rbind(results, temp_results)
    }
    nComps <- min(results[results$P.value > 0.05, "Component"])
    print(paste0("*** Optimal number of components based on t.test: ", nComps))
    bp <- ggplot(pressDFres, aes(x = variable, y = value)) + 
      theme_bw() + geom_boxplot(notch = FALSE) + labs(x = "Number of Components", 
                                                      y = "PRESS") + 
      stat_boxplot(geom = "errorbar", width = 0.2) + 
      geom_vline(xintercept = nComps, linetype = "dashed", 
                 color = "blue", size = 1)
    theme(axis.text = element_text(size = 18), legend.position = "none", 
          axis.title = element_text(size = 20, face = "bold"), 
          axis.text.x = element_text(angle = 0, vjust = 0.5), 
          panel.border = element_rect(linetype = "solid", 
                                      fill = NA, size = 1.5))
    print(bp)
  }
  if (method == "firstMin") {
    press.out <- spectratrait::pls_permutation_by_groups(dataset = dataset, targetVariable=inVar,
                                                   maxComps=maxComps, iterations=iterations, 
                                                   prop=prop, group_variables=group_variables)
    pressDF <- as.data.frame(press.out$PRESS)
    names(pressDF) <- as.character(seq(maxComps))
    pressDFres <- reshape2::melt(pressDF)
    mean_PRESS_comp <- apply(X = pressDF, MARGIN = 2, FUN = mean)
    lowest_PRESS <- which.min(mean_PRESS_comp)
    results <- as.vector(array(data = "NA", dim = c(lowest_PRESS - 1, 1)))
    for (i in seq_along(1:(lowest_PRESS - 1))) {
      comp1 <- i
      comp2 <- lowest_PRESS
      ttest <- t.test(pressDFres$value[which(pressDFres$variable == comp1)], 
                      pressDFres$value[which(pressDFres$variable == comp2)])
      results[i] <- round(unlist(ttest$p.value), 8)
    }
    results <- data.frame(seq(1, lowest_PRESS - 1, 1), results)
    names(results) <- c("Component", "P.value")
    first <- min(which(as.numeric(as.character(results$P.value)) > 0.05))
    nComps <- results$Component[first]
    print(paste0("*** Optimal number of components based on t.test: ", nComps))
    bp <- ggplot(pressDFres, aes(x = variable, y = value)) + 
      theme_bw() + geom_boxplot(notch = FALSE) + labs(x = "Number of Components", 
                                                      y = "PRESS") + 
      stat_boxplot(geom = "errorbar", width = 0.2) + 
      geom_vline(xintercept = nComps, linetype = "dashed", 
                 color = "blue", size = 1)
    theme(axis.text = element_text(size = 18), legend.position = "none", 
          axis.title = element_text(size = 20, face = "bold"), 
          axis.text.x = element_text(angle = 0, vjust = 0.5), 
          panel.border = element_rect(linetype = "solid", 
                                      fill = NA, size = 1.5))
    print(bp)
  }
  return(nComps)
}