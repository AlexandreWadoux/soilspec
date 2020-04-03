\name{mineralRef}
\alias{mineralRef}
\title{Reference spectra for clay minerals and iron oxides}
\usage{
data('mineralRef')
}
\description{
Specific mineral compounds measured by laboratory spectrometers, from Kokaly et al., (2017).
A data frame with 2151 rows and 13 columns:
\itemize{
\item{\bold{wavelength}: wavelength (in nm) of the reference spectra.}
\item{\bold{goethite}: reference spectrum for goethite.}
\item{\bold{hematite}: reference spectrum for hematite. }
\item{\bold{gypsum}: reference spectrum for gypsum.}
\item{\bold{calcite}: reference spectrum for calcite.}
\item{\bold{kaolinite}: reference spectra for kaolinite. Two spectra are provided: \code{kaolinite_114} and \code{kaolinite_113}.}
\item{\bold{montmorillonite}: reference spectra for montmorillonite. Two spectra are provided: \code{montmorillinite_126} and \code{montmorillinite_127}.}
\item{\bold{illite}: reference spectra for illite. Two spectra are provided: \code{illite_121} and \code{illite_120}.}
\item{\bold{kaolinite/smectite}: : reference spectra for the mixture kaolinite/smectite. Two spectra are provided: \code{kaol_smect_124} and \code{kaol_smect_125}.}
}
}
\details{
The reference spectra are freely available by the USGS through their website. The minerals are purified material so that each reference spectra contain unique spectro-chemical links related to chemical structure of the mineral. The collected reference spectra are for clay minerals including kaolinite (KGa-2), illite (GDS4),smectite (SWy-1), kaolinite-smectite 50/50 mixture (H89-FR-2), goethite (GDS240), and hematite (GDS576).
}

\examples{
## load the data into the environment
data('mineralRef')

## take the reference spectra for calcite
calcite <- mineralRef$calcite
}
\references{
Kokaly, R.F., Clark, R.N., Swayze, G.A., Livo, K.E., Hoefen, T.M., Pearson, N.C., ... & Klein, A.J. (2017). USGS spectral library version 7 (No. 1035). US Geological Survey.
}
