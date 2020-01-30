\name{chBLC_ext}
\alias{chBLC_ext}
\title{Convex hull function used for regions of interest}
\usage{
chBLC_ext(spectra, lower, upper, ...)
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
  \item{c.hull}{the fitted convex hull to the region of interest.}
  \item{raw.spec}{the raw spectra used for trimming.}
  \item{continuum}{the continuum removed spectra.}
  \item{polygon}{the polygon of the area between the continnum removed spectra and the line of the convex hull.}
  }

\examples{

require(soilspec)

# load the reference spectra
data("mineral_ref")

# inspect the object
str(mineral_ref)

# create a sequence of number to represent the wavelength
wavelength <- seq(350, 2500, by = 1)

#Function for trimming spectra or targeting a specific spectral region of #interest

trimSpec <- function(spectra, wavlimits) {
  datawavs <- as.numeric(names(spectra))
  limits <- which(datawavs \%in\% wavlimits)
  kept_index <- seq(limits[1], limits[2], 1)
  trimmed_spectra <- spectra[, kept_index]
  kept_names <- datawavs[kept_index]
  colnames(trimmed_spectra) <- kept_names
  return(trimmed_spectra)}

# select reference mineral
mont_ref <- as.data.frame(t(mineral_ref$kaolinite_114))
colnames(mont_ref) <- wavelength

# diagnostic wavelengths
lower <- 2078
upper <- 2267

# spectrum trimming
mont_diog1 <- trimSpec(mont_ref,
                       wavlimits = range(lower:upper))
mont_diog1.CR <- chBLC_ext(spectra = mont_diog1, lower = lower, upper = upper)

# plot the region of interest and the fitted convex hull to the continuum removed spectra
plot(seq(from = lower, to = upper, by = 1), mont_diog1[1,],
     type = "l",
     xlab = "Wavelength /nm",
     ylab = "Reflectance")

# continuum removed spectra
lines(seq(from = lower, to = upper, by = 1), mont_diog1.CR$continuum,
      col = "red")


feature.bd <- 1 - (mont_diog1.CR[[3]][which(mont_diog1.CR[[2]]==min(mont_diog1.CR[[2]]))]/mont_diog1.CR[[4]][which(mont_diog1.CR[[2]]==min(mont_diog1.CR[[2]]))])

# display band depth
feature.bd

}
