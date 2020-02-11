\name{chBLCext}
\alias{chBLCext}
\title{Convex hull function used for regions of interest}
\usage{
chBLCext(spectra, lower, upper, ...)
}
\arguments{
  \item{spectra}{Pre-processed spectra (not trimmed).}
  \item{lower}{Lower limit of the region of interest.}
  \item{upper}{Upper limit of the region of interest.}
  \item{...}{Additional arguments to be passed to the function.}
}
\description{
This function fits a convex hull to the region of interest. The regions of interest are absorbance (or reflectance) of some secondary clay minerals or iron oxides from vis-NIR spectra.}
\references{
Malone, B. P., Hughes, P., McBratney, A. B., & Minasny, B. (2014). A model for the identification of terrons in the Lower Hunter Valley, Australia. Geoderma Regional, 1, 31-47.
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
montDiog1CR <- chBLCExt(spectra = montDiog1, lower = lower, upper = upper)

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
