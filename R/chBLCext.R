#' @title chBLCext
#'
#' @description Fitting a Convex hull function for regions of interest
#'
#' @param spectra A matrix of spectra
#'
#' @return the fitted convex hull
#'
#' @export

chBLCext <- function(spectra, lower, upper, type, ...){
  
  wav <- as.numeric(colnames(spectra))
  hullSpectra <- matrix(NA, ncol = ncol(spectra), nrow = nrow(spectra))
   tempSpect = as.matrix(spectra)
  # data1 <- sortedXyData(wav, tempSpect)
  
  X <- spectra
  cHullfun <- function(x, wav) {
    id <- sort(chull(c(wav[1] - 1, wav, wav[length(wav)] + 1), c(0, x, 0)))
    cHull <- id[-c(1, length(id))] - 1
    return(cHull)
  }
  
  if (type == "A") {
    X <- 1/X
}

cHull <- cHullfun(X, wav)

## calculate linear approximation between hull points
linearApprox <- approx(x = wav[cHull], y = tempSpect[cHull] , xout = wav, method = 'linear', ties = 'mean')
## calculate the deviation from the convex hull
hullSpectra[1,] <- 1 - (( linearApprox[[2]] - tempSpect )/linearApprox[[2]])
colnames(hullSpectra) <- wav
# construct polygon
aaa <- hullSpectra
aaa[,c(1,ncol(aaa))] <- 1.0001
bbb <- cbind(wav, t(aaa))
ccc <- cbind(wav, rep(c(1.0001), each = length(wav)))
ddd <- as.matrix(rbind(bbb, ccc))
retval <- list(wave = wav, cHull = hullSpectra, rawSpec = tempSpect, continuum = linearApprox[[2]], polygon = ddd)
return(retval)

}
