\name{spectra2colour}
\alias{spectra2colour}
\title{Function to extract RGB and Munsell colour from a spectra}
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
A dataframe with of size the number of bands plus two. The columns correspond to the number of spectral bands. The RBG colours are accessed by the column name `colour` and the Munsel colours are accessed by the column name `munsell`. One colour is provided for each spectra. The spectra is also returned.
}
