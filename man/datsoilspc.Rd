\name{datsoilspc}
\alias{datsoilspc}
\title{Geeves dataset}
\usage{
data('datsoilspc')
}
\description{

The Geeves dataset (Geeves et al., 1994) is a soil survey containing 391 soil samples collected in the wheatbelt of southern New South Wales and northern Victoria, Australia. The reflectance infrared spectra was obtained in the range 350-2500 nm at 1 nm resolution using an AgriSpec instrument with a contact probe (Analytical Spectral Devices, Boulder, Colorado, USA). A spectralon (Labsphere Inc., North Sutton, N.H., USA) was used as reflectance standard. Each spectrum in the dataset is a wavelength-average of 40 scans.

Total carbon (TC), clay, silt and sand were analysed in the laboratory using standard methods. TC was analysed by the dry-combustion method and is reported in g 100 g-1. The values of the TC range from 0.6 to 12.74 g 100 g-1. The TC does not contain inorganic carbon. The hydrometer method was used for the particle size distribution and their relative quantity was converted to percent. Clay content varies from 5 to 74\%, silt content from 2 to 54\% and sand content from 14 to 91\%.
\itemize{
   \item{\bold{clay}: Clay content (\%).}
   \item{\bold{silt}: Silt content (\%).}
   \item{\bold{sand}: Sand content (\%).}
   \item{\bold{TotalCarbon}: Total carbon g 100 g-1.}
   \item{\bold{spc}: Table of size 391 by 2191 with vis-NIR spectra. Each row is a spectra.}
}
}
\references{
Geeves G.W., Cresswell H.P., Murphy B.W., Gessler P.I., Chartres C.J., Little I.P., Bowman G.M. (1994) Physical, chemical and morphological properties of soils in the wheat-belt of southern NSW and northern Victoria. NSW Department of Conservation; Land Management/CSIRO Division of Soils Occasional Report, CSIRO, Australia
}
\examples{
## load the data into the environment
data('datsoilspc')
}
