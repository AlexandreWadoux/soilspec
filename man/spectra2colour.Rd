\name{spectra2colour}
\alias{spectra2colour}
\title{Function to extract RGB and Munsell color from a spectra}
\usage{
spectra2colour(spectra, ...)
}
\arguments{
\item{spectra}{a matrix or dataframe containing the spectra}
\item{...}{other arguments to be passed in the function}
}
\description{
Converts spectra reflectance into RGB and Munsell colours
}
\value{
A dataframe with of size the number of bands plus two. The columns correspond to the number of spectral bands. The RBG colors are accessed by the column name `colour` and the Munsel color are accessed by the column name `munsell`. One color is provided for each spectra. The spectra is also returned.
}
