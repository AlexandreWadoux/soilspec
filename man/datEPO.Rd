\name{datEPO}
\alias{datEPO}
\docType{data}
\title{Spectra with moisture}
\description{
Subset of the Geeves soil survey data in which the soil samples were scanned at different level of water content.  
\code{datEPO} is a list with four objects:
\itemize{
   \item \code{soilC}: benchmark spectra, data frame with 100 rows and 1 column.
   \item \code{spectra0}: 100 absorbance spectra from air-dried soil, data frame with 100 rows and 2151 columns.
   \item \code{spectra1}: 100 absorbance spectra from wet-soil, data frame with 100 rows and 2151 columns.
   \item \code{spectra2}: 100 absorbance spectra from wet-soil after being air-dried for 1 day, data frame with 100 rows and 2151 columns.
}
}
\usage{data("datEPO")}

\references{
Minasny, B. et al. (2011) Removing the effect of soil moisture from NIR diffuse reflectance spectra for the prediction of soil organic carbon. Geoderma, 167, 118-124.}
\examples{
data(datEPO)

soilC <- datEPO$soilC
spectra0 <- datEPO$spectra0
spectra1 <- datEPO$spectra1
spectra2 <- datEPO$spectra2

}
\keyword{datasets}
