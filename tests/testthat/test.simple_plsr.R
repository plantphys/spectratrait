context("*** Test simple PLSR model fit using Ely et al. (2019) data ***")

### Setup data for tests
#Load Ely et al 2019 dataset
data("ely_plsr_data")
inVar <- "N_g_m2"
Start.wave <- 500
End.wave <- 2400
wv <- seq(Start.wave,End.wave,1)
plsr_data <- ely_plsr_data
spec <- as.matrix(plsr_data[, which(names(plsr_data) %in% paste0("Wave_",wv))])
plsr_data <- data.frame(plsr_data[, which(names(plsr_data) %notin% paste0("Wave_",wv))],
                        Spectra=I(spec))

test_that("Test fitting simple PLSR model", {
  method <- "firstPlateau"
  random_seed <- 1245565
  seg <- 50
  maxComps <- 20
  iterations <- 80
  prop <- 0.70
  
  nComps <- spectratrait::find_optimal_components(dataset=plsr_data, targetVariable=inVar, 
                                                  method=method, maxComps=maxComps, 
                                                  iterations=iterations, seg=seg, prop=prop, 
                                                  random_seed=random_seed)
  
  segs <- 100
  plsr.out <- plsr(as.formula(paste(inVar,"~","Spectra")),scale=FALSE,ncomp=nComps,validation="CV",
                   segments=segs, segment.type="interleaved",trace=FALSE,data=plsr_data)
  fit <- plsr.out$fitted.values[,1,nComps]
  expect_type(fit, "double")
})