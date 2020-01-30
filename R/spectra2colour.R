#' Method that extracts RGB and Munsell colour from a matri or dataframe of spectra.
#'
#' Converts spectra reflectance into RGB and Munsell colours
#' @importFrom plyr adply
#' @importFrom plyr splat
#' @importFrom plyr llply
#' @importFrom munsell rgb2mnsl
#'
#' @param spectra a matrix or dataframe containing the spectra
#' @param ... other arguments to be passed in the function
#'
#' @exportMethod spectra2colour
#' @author Michael Nelson and Mario Fajardo
spectra2colour <- function(spectra, ...){

  #Internal Functions
  spectra_to_RGB <- function(.spectra,
                             all_wavelengths,
                             rgb_list = list(blue = list(.interval=450:520),
                                             red = list(.interval=600:690),
                                             green=list(.interval=520:600)),
                             ...) {
    ## a function to return the average values in the
    ## red, green and blue sections of the spectrum
    ## would work on any intervals
    ##
    ## get the appropriate indices
    interval_list <- plyr::llply(rgb_list, plyr::splat(in_interval), .all=all_wavelengths)
    ##
    ## get the average in these subsets
    rgb_values <- lapply(interval_list, mean_interval, .data=.spectra)
    ##
    ## convert to colour
    colour <- with(rgb_values, rgb(red, green, blue))
    ##
    ## return data frame
    with(rgb_values, data.frame(red, green, blue, colour))
  }

  in_interval <- function(.all, .interval,...){
    ## index for an interval
    ## a function to subset a particular waveband interval
    .in_interval = .all %in% .interval
  }

  mean_interval <- function(.data, .index){
    ## returns the mean for given indices
    mean(as.numeric(.data[.index]), na.omit = TRUE)
  }

  ##End of internal functions


  ## find r,g,b colour
  rgb_colours <- plyr::adply(spectra, 1, spectra_to_RGB, all_wavelengths = colnames(spectra),.id = NULL)
  ##
  ## get munsell colour
  munsell_colours <- plyr::splat(function(red,green,blue, ...){munsell::rgb2mnsl(R=red,G=green,B=blue)})(rgb_colours)
  ##
  ## return
  soil_colours <- data.frame(rgb_colours, munsell = munsell_colours)

  return(soil_colours)
}

