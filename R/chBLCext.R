#' @title chBLCext
#'
#' @description Fitting a Convex hull for regions of interest
#'
#' @param spectra a matrix containing at least one spectrum of a given region 
#' of interest
#' @param type FIXME!
#' @param wav (optional) numeric vector containing the wavelengths of the region 
#' @return a list with the following components 
#' \itemize{
#'   \item{\code{\link{spectra}}: the input \code{spectra}}
#'   \item{\code{\link{wav}}: the input \code{wav}}
#'   \item{\code{\link{continuum_removed}}: a matrix with the removed continuum}
#'   \item{\code{\link{continuum}}: a matrix with the continuum spectra}
#'   \item{\code{\link{polygon}}: FIXME! not sure what this is}
#'}
#' @export
#------------------------------- Info ------------------------------------------
# Function inspired by the 'continuumRemoval' function of the prospectr package
# Author: Alexandre Wadoux
#
# Date:        April 2020
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
chBLCext <- function(spectra, type = c("R", "A"), wav) {

  # ... (ellipsis) is not necessary in this function as no unspecified arguments
  # being passed to any internal function

  # NOTE that the simple inverse might cause problems in some spectra,
  # instead transform to Reflectance
  if (type == "A") {
    spectra <- 1 / 10^spectra
  }
  spectra <- cbind(0, spectra, 0)

  # make wav optional
  if (missing(wav)) {
    wav <- 1:ncol(spectra)
    output_wav <- NULL
  } else {
    output_wav <- wav
  }

  wav <- c(wav[1] - 1, wav, wav[length(wav)] + 1)

  # define function get_chull_spc() to get the hull spectrum of a single
  # input spectrum/vector
  get_chull_spc <- function(x, y = wav) {
    hull_spc <- sort(chull(x, y))
    hull_spc[1] <- hull_spc[1] + 1
    hull_spc[length(hull_spc)] <- hull_spc[length(hull_spc)] - 1
    hull_spc <- unique(hull_spc)
    hull_interp <- approx(
      x = y[hull_spc],
      y = x[hull_spc],
      xout = y,
      method = "linear"
    )$y
    return(hull_interp)
  }

  spc_hull <- t(apply(spectra, MARGIN = 1, FUN = get_chull_spc))
  spc_hull <- spc_hull[, -c(1, ncol(spc_hull))]

  spectra <- spectra[, -c(1, ncol(spectra))]
  wav <- wav[-c(1, length(wav))]

  if (type == "A") {
    # substraction - absorbance (Fig. 5 Clark & Roush (1984))
    hullSpectra <- 1 + spectra - spc_hull
  } else {
    # division - reflectance (Fig. 5 Clark & Roush (1984))
    hullSpectra <- spectra / spc_hull
  }

  if (type == "A") {
    # back transform the CR spectra
    hullSpectra <- log(1 / hullSpectra, 10)
    # back transform to show the CH line
    spc_hull <- log(1 / spc_hull, 10)
    # back transform the spectra
    spectra <- log(1 / spectra, 10)
  }

  # prepare xy for the polygon
  pol <- cbind(wav, as.numeric(hullSpectra[1, ]))

  retval <- list(
    spectra = spectra,
    wav = output_wav,
    continuum_removed = hullSpectra,
    continuum = spc_hull,
    polygon = pol
  )
  return(retval)
}
