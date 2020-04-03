
\name{datEPO}
\alias{datEPO}
\docType{data}
\title{Spectra scanned on moist soil samples}
\description{
Subset of the Geeves soil survey data in which the soil samples were scanned with different level of water content. List with four data frames.

\itemize{
   \item{\bold{soilC}: total carbon values, in g 100 g-1, data frame with 100 rows and 1 column.}
   \item{\bold{spectra0}: 100 absorbance spectra from air-dried soil, data frame with 100 rows and 2151 columns (average moisture 5\%).}
   \item{\bold{spectra1}: 100 absorbance spectra from wet-soil, data frame with 100 rows and 2151 columns (average moisture 12\%).}
   \item{\bold{spectra2}: 100 absorbance spectra from wet-soil after being air-dried for one day, data frame with 100 rows and 2151 columns (average moisture 9\%).}
}
}
\usage{data("datEPO")}

\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Minasny, B., McBratney, A.B., Bellon-Maurel, V., Roger, J.M., Gobrecht, A., Ferrand, L., & Joalland, S. (2011). Removing the effect of soil moisture from NIR diffuse reflectance spectra for the prediction of soil organic carbon. Geoderma, 167, 118-124.}
\examples{
## load the dataset
data(datEPO)

## take object from the list
soilC <- datEPO$soilC
spectra0 <- datEPO$spectra0
spectra1 <- datEPO$spectra1
spectra2 <- datEPO$spectra2
