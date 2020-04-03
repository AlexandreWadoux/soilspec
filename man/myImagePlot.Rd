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
Minasny, B., McBratney, A.B., Bellon-Maurel, V., Roger, J.M., Gobrecht, A., Ferrand, L., & Joalland, S. (2011). Removing the effect of soil moisture from NIR diffuse reflectance spectra for the prediction of soil organic carbon. Geoderma, 167, 118-124.
}
