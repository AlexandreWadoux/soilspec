\name{datsoilspc}
\alias{datsoilspc}
\title{Geeves dataset}
\usage{
data('datsoilspc')
}
\description{

The Geeves dataset (Geeves et al., 1994) is a soil survey containing 391 soil samples collected in the wheatbelt of southern New South Wales and northern Victoria, Australia. The reflectance infrared spectra was obtained in the range 350-2500 nm at 1 nm resolution using an AgriSpec instrument with a contact probe (Analytical Spectral Devices, Boulder, Colorado, USA). A spectralon (Labsphere Inc., North Sutton, N.H., USA) was used as reflectance standard. Each spectrum in the dataset is a wavelength-average of 40 scans.

Total carbon (TC), clay, silt and sand were analysed in the laboratory using standard methods. TC was analysed by the dry-combustion method and is reported in $\mathrm{g~100~g}^{-1}$. The values of the TC range from 0.6 to 12.74 \eqn{\mathrm{g~100~g}^{-1}}. The TC does not contain inorganic carbon. The hydrometer method was used for the particle size distribution and their relative quantity was converted to percent. Clay content varies from 5 to 74%, silt content from 2 to 54% and sand content from 14 to 91%.
   
   \item{clay}{Clay content (%)}
   \item{silt}{Silt content (%)}
   \item{sand}{Sand content (%)}
   \item{TotalCarbon}{Total carbon (\eqn{\mathrm{g~100~g}^{-1}}.)}
   \item{spc}{Table of size 391 \eqn{\times} 2191 with vis-NIR spectra. Each row is a spectra.}   
}
\references{
Geeves GW, Cresswell HP, Murphy BW, Gessler PI, Chartres CJ, Little IP, Bowman GM (1994) Physical, chemical and morphological properties of soils in the wheat-belt of southern NSW and northern Victoria. NSW Department of Conservation; Land Management/CSIRO Division of Soils Occasional Report, CSIRO, Australia
}
\examples{
data('datsoilspc')
}
