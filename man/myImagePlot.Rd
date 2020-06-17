\name{myImagePlot}
\alias{myImagePlot}
\title{Plot a projection matrix}
\usage{
myImagePlot(x, ...)
}
\arguments{
  \item{x}{A projection matrix}

  \item{...}{Other arguments to be passed in the function}
}
\description{
Function to plot an image from a projection matrix.
}
\examples{
## create a diagonal matrix
diagM <- diag(1, 100, 100)

## plot the diagonal matrix
myImagePlot(diagM)

}
\value{
Return a plot the input matrix.
}
\references{
Mehmood, T., Liland, K. H., Snipen, L., & Saebo, S. (2012). A review of variable selection methods in partial least squares regression. Chemometrics and Intelligent Laboratory Systems, 118, 62-69.
