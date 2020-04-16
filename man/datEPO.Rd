\name{datEPO}
\alias{datEPO}
\docType{data}
\title{Spectra with moisture}
\arguments{
A list called `datEPO` with four objects:
   \item{soilC}{benchmark spectra, data frame with 100 rows and 1 column.}
   \item{spectra0}{100 absorbance spectra from air-dried soil, data frame with 100 rows and 2151 columns.}
   \item{spectra1}{100 absorbance spectra from wet-soil, data frame with 100 rows and 2151 columns.}
   \item{spectra2}{100 absorbance spectra from wet-soil after being air-dried for 1 day, data frame with 100 rows and 2151 columns.}

}
\usage{data("datEPO")}

\details{
Subset of the Geeves soil survey data in which the soil samples were scanned at different level of water content. Data for the Section 'Removing soil moisture (EPO)'

}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
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
