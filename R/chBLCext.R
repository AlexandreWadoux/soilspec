#' @title chBLCext
#'
#' @description Fitting a Convex hull function for regions of interest
#'
#' @param spectra A matrix of spectra
#'
#' @return the fitted convex hull
#'
#' @export

#Convex hull function used for regions of interest work such as detecting secondary clay minerals from vis-NIR spectra
chBLCext <- function(spectra, lower, upper, ...){
  interval <- c(1:ncol(spectra))
  hullSpectra <- matrix(NA, ncol = ncol(spectra), nrow = nrow(spectra))
  tempSpect = as.matrix(spectra)
  data1 <- sortedXyData(interval, tempSpect)

  ## calculate convex hull
  cHull <- chull(data1)
  ## get the appropriate region
  cc <-cHull
  xs <- which(cc == 1)
  ccd <- cc[which(cc == 1):length(cc)]
  ccdM <- matrix(NA, nrow = 1, ncol = length(ccd))
  for(f in 1:(length(ccd)-1)){
    if(ccd[f] < ccd[f+1]){ccdM[1,f] <- ccd[f]}
    else {ccdM[1,f] <- ccd[f]
    break}}
  if(f == (length(ccd)-1)){ccdM[1,ncol(ccdM)] <- ccd[length(ccd)]}
  tempX <- c(ccdM)[!is.na(c(ccdM))]
  #Left side
  if(cHull[1] > cHull[length(cHull)]){lss <- cc[1:(xs - 1)]
  slss <- lss[which(lss >= lss[1])]
  cHull <- sort(c(tempX, slss))} else {cHull <- tempX}

  ## calculate linear approximation between hull points
  linearApprox <- approx(data1[cHull,], xout = interval, method = 'linear', ties = 'mean')
  ## calculate the deviation from the convex hull
  hullSpectra[1,] <- 1- (( linearApprox[[2]] - tempSpect )/linearApprox[[2]])
  colnames(hullSpectra) <- colnames(spectra)
  # construct polygon
  aaa <- hullSpectra
  aaa[,c(1,ncol(aaa))] <- 1.0001
  bbb <- cbind(seq(lower, upper, by = 1), t(aaa))
  ccc <- cbind(seq(lower, upper,by = 1), rep(c(1.0001),each = length(seq(lower, upper, by = 1))))
  ddd <- as.matrix(rbind(bbb, ccc))
  retval <- list(wave = seq(lower, upper, by = 1),cHull = hullSpectra, rawSpec = data1$y,continuum = linearApprox[[2]], polygon = ddd)
  return(retval)}
