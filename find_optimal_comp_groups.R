##Same function as find_optimal_components but just with an additional argument (groups)
## and using the new function pls_permut_by_group instead of pls_permutation

##Stratifying is only included in firstMin and firstPlateau, maybe pls could be eliminated...
##or otherwise some way of stratifying should be included in the pls method (I don't know how).


## groups should be an vector of caracters of the form c("var1", "var2"..."varn")
##I have not included any control for this condition, though

##Changed commands marked 

find_optimal_comp_groups <- function (dataset = NULL, method = "pls", maxComps = 20, 
          iterations = 20, seg = 100, prop = 0.7, random_seed = 123456789, groups=NULL) 
{
  set.seed(random_seed)
  if (method == "pls") {
    print("*** Running PLS permutation test ***")
    plsr.out <- pls::plsr(as.formula(paste(inVar, "~", 
                                           "Spectra")), scale = FALSE, center = TRUE, 
                          ncomp = maxComps, validation = "CV", segments = seg, 
                          segment.type = "interleaved", trace = FALSE, 
                          jackknife = TRUE, data = cal.plsr.data)
    nComps <- selectNcomp(plsr.out, method = "onesigma", 
                          plot = TRUE)
  }
  if (method == "firstPlateau") {
##Changed pls_permutation for pls_permut_by_groups    
    press.out <- spectratrait::pls_permut_by_groups(dataset = dataset, 
                                               maxComps = maxComps, iterations = iterations, prop = prop, groups=groups)
    pressDF <- as.data.frame(press.out$PRESS)
    names(pressDF) <- as.character(seq(maxComps))
    pressDFres <- reshape2::melt(pressDF)
    results <- NULL
    for (i in 1:(maxComps - 1)) {
      p_value <- t.test(press.out$PRESS[, i], press.out$PRESS[, 
                                                              (i + 1)])$p.value
      temp_results <- data.frame(Component = (i + 1), P.value = round(p_value, 
                                                                      6))
      results <- rbind(results, temp_results)
    }
    nComps <- min(results[results$P.value > 0.05, "Component"])
    print(paste0("*** Optimal number of components based on t.test: ", 
                 nComps))
    bp <- ggplot(pressDFres, aes(x = variable, y = value)) + 
      theme_bw() + geom_boxplot(notch = FALSE) + labs(x = "Number of Components", 
                                                      y = "PRESS") + stat_boxplot(geom = "errorbar", 
                                                                                  width = 0.2) + geom_vline(xintercept = nComps, linetype = "dashed", 
                                                                                                            color = "blue", size = 1)
    theme(axis.text = element_text(size = 18), legend.position = "none", 
          axis.title = element_text(size = 20, face = "bold"), 
          axis.text.x = element_text(angle = 0, vjust = 0.5), 
          panel.border = element_rect(linetype = "solid", 
                                      fill = NA, size = 1.5))
    print(bp)
  }
  if (method == "firstMin") {
##Changed pls_permutation for pls_permut_by_groups
    press.out <- spectratrait::pls_permut_by_goups(dataset = dataset, 
                                               maxComps = maxComps, iterations = iterations, prop = prop, groups=groups)
    pressDF <- as.data.frame(press.out$PRESS)
    names(pressDF) <- as.character(seq(maxComps))
    pressDFres <- reshape2::melt(pressDF)
    mean_PRESS_comp <- apply(X = pressDF, MARGIN = 2, FUN = mean)
    lowest_PRESS <- which.min(mean_PRESS_comp)
    results <- as.vector(array(data = "NA", dim = c(lowest_PRESS - 
                                                      1, 1)))
    for (i in seq_along(1:(lowest_PRESS - 1))) {
      comp1 <- i
      comp2 <- lowest_PRESS
      ttest <- t.test(pressDFres$value[which(pressDFres$variable == 
                                               comp1)], pressDFres$value[which(pressDFres$variable == 
                                                                                 comp2)])
      results[i] <- round(unlist(ttest$p.value), 8)
    }
    results <- data.frame(seq(1, lowest_PRESS - 1, 1), results)
    names(results) <- c("Component", "P.value")
    first <- min(which(as.numeric(as.character(results$P.value)) > 
                         0.05))
    nComps <- results$Component[first]
    print(paste0("*** Optimal number of components based on t.test: ", 
                 nComps))
    bp <- ggplot(pressDFres, aes(x = variable, y = value)) + 
      theme_bw() + geom_boxplot(notch = FALSE) + labs(x = "Number of Components", 
                                                      y = "PRESS") + stat_boxplot(geom = "errorbar", 
                                                                                  width = 0.2) + geom_vline(xintercept = nComps, linetype = "dashed", 
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