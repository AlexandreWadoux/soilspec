#' @title chBLCext
#'
#' @description Fitting a Convex hull for regions of interest
#'
#' @param a spectrum 
#'
#' @return the fitted convex hull
#'
#' @export
#------------------------------- Info ------------------------------------------
# Function inspired by the 'continuumRemoval' function of the prospectr package
# Author: Alexandre Wadoux            
#
# Date:        April 2020
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
chBLCext<- function(spectra, type = c("R", "A"), wav, ...){
  
  X <- spectra
  
  # if absorbance instead take the inverse to find the CH points
  if (type == "A"){X <- 1/X}
  
  cHullFun <- function(x, wav, interpol) {
    cHull <- sort(chull(c(wav[1] - 1, wav, wav[length(wav)] + 1), c(0, x, 0)))
    cHull <- cHull[-c(1, length(cHull))] - 1
    cont <- approx(x = wav[cHull], y = x[cHull], xout = wav, method = "linear")$y
    return(cont)
  }
  
  cont <- cHullFun(X, wav, interpol)
  
  if (type == "A"){ 
    # substraction - absorbance (Fig. 5 Clark & Roush (1984))
    hullSpectra <-  1 + X - cont 
  }else{
    # division - reflectance (Fig. 5 Clark & Roush (1984))
    hullSpectra <-  X/cont  
  }
  
  if (type == "A"){ 
    # back transform the CR spectra
    hullSpectra <- 1/hullSpectra  
    # back transform to show the CH line
    cont <- 1/cont 
    # back transformt the spectra
    X <- 1/X 
  }
  
  # prepare xy for the polygon
  pol <- cbind(wav, as.numeric(hullSpectra[1,]))
  
  
  retval<- list(wave = wav, 
                cHull = hullSpectra, 
                rawSpec = X, 
                continuum = cont, 
                polygon = pol)
  return(retval)
}

