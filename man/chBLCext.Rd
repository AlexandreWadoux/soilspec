\name{chBLCext}
\alias{chBLCext}
\title{Convex hull function used for regions of interest}
\usage{
chBLCext(spectra, type, wav, ...)
}
\arguments{
  \item{spectra}{Pre-processed spectrum for the region of interest. Only one spectrum can be provided.}
  \item{type}{Either 'A' of the spectrum is in absorbance, or 'R' if the spectrum is in reflectance.}
  \item{wav}{Vector numeric containing the wavelengths the of \code{spectra} argument.}

  \item{...}{Additional arguments to be passed to the function.}
}
\description{
This function fits a convex hull to the region of interest. The regions of interest are reflectance of absorbance of some secondary clay minerals or iron oxides from vis-NIR spectra. Only one spectrum can be provided, not a spectra table.}
\references{
Clark, R.N., and Roush, T.L., 1984. Reflectance Spectroscopy: Quantitative Analysis Techniques for Remote Sensing Applications. J. Geophys. Res. 89, 6329-6340.
  }
\value{
A list with the following arguments:
  \item{wave}{a vector of the wavelengths for the trimmed spectra.}
  \item{cHull}{the fitted convex hull to the region of interest.}
  \item{rawSpec}{the raw spectra used for trimming.}
  \item{continuum}{the continuum removed spectra.}
  \item{polygon}{the polygon of the area between the continnum removed spectra and the line of the convex hull.}
  }

\examples{

require(soilspec)

# load the reference spectra
data("mineralRef")

# inspect the object
str(mineralRef)

# create a sequence of number to represent the wavelength
wavelength <- seq(350, 2500, by = 1)

#Function for trimming spectra or targeting a specific spectral region of #interest

trimSpec <- function(spectra, wavlimits){
  datawavs <- as.numeric(colnames(spectra))
  limits <- which(datawavs %in% wavlimits)
  keptIndex <- seq(limits[1], limits[2], 1)
  trimmedSpectra <- spectra[, keptIndex]
  keptNames <- datawavs[keptIndex]
  colnames(trimmedSpectra) <- keptNames
  return(trimmedSpectra)
}

# select reference mineral
montRef <- as.data.frame(t(mineralRef$kaolinite_114))
colnames(montRef) <- wavelength

# diagnostic wavelengths
lower <- 2078
upper <- 2267

# spectrum trimming
montDiog1 <- trimSpec(montRef,
                       wavlimits = range(lower:upper))
montDiog1CR <- chBLCext(spectra = montDiog1, wav = wavelength, type = 'R')

# plot the region of interest and the fitted convex hull to the continuum removed spectra
plot(seq(from = lower, to = upper, by = 1), montDiog1[1,],
     type = "l",
     xlab = "Wavelength /nm",
     ylab = "Reflectance")

# continuum removed spectra
lines(seq(from = lower, to = upper, by = 1), montDiog1CR$continuum,
      col = "red")

# band depth
featureBd <- 1 - (montDiog1CR[[3]][which(montDiog1CR[[2]]==min(montDiog1CR[[2]]))]/montDiog1CR[[4]][which(montDiog1CR[[2]]==min(montDiog1CR[[2]]))])
featureBd

}
