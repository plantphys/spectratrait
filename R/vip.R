##' @title VIP returns all VIP values for all variables and all number of components, as a ncomp x nvars matrix.
##' @param object fitted pls::plsr object
##' @export
VIP <- function(object) {
  ## VIP returns all VIP values for all variables and all number of components,
  ## as a ncomp x nvars matrix.
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*") # Replace with matrix mult.
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}

##' @title VIPjh returns the VIP of variable j with h components
##' @param object fitted pls::plsr object
##' @param j which variable in the fitted pls::plsr object
##' @param h the number of components in the fitted pls::plsr object to calculate the VIP
##' @export
VIPjh <- function(object, j, h) {
  ## VIPjh returns the VIP of variable j with h components
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