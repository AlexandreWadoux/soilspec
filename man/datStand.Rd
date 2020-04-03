\name{datStand}
\alias{datStand}
\title{Vis-NIR spectra of a standard material}
\usage{
data('datStand')
}
\description{
Vis-NIR spectra of a standard material as described in Ben-Dor et al., (2015). The specific sample is the "Lucky Bay" sample that is essentially bright (white colored) sand material collected from a homogeneous sand dune situated at Lucky Bay in southwestern Australia. This standard material can be requested from the authors of Ben-Dor et al., (2015) which comes with an associated vis-NIR spectrum measured with their instrumentation and to their published protocol, and used to benchmark other spectrometers.

\itemize{
   \item{\bold{CSIRO.s}: benchmark LB spectra, data frame with two rows and 2151 columns.}
   \item{\bold{ISS_slaves}: spectra collected from two instruments, data frame with two rows and 2151 columns.}
   \item{\bold{asd.soil}: two soil samples scans by an ASD instrument, data frame with two rows and 2151 columns.}
   \item{\bold{se.soil}: two soil samples scans by an SE instrument, data frame with two rows and 2151 columns.}
}
}

\references{
Ben-Dor, E., Ong, C., & Lau, I.C. (2015). Reflectance measurements of soils in the laboratory: Standards and protocols. Geoderma, 245, 112-124.
}
\examples{
## load the data into the environment
data('datStand')
}
