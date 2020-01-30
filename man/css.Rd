\name{css}
\alias{css}
\title{Adequate calibration set size for infrared spectroscopy data}
\usage{
css(S, k, method = "kms", repetitions = 10, n = 512, from, to, bw, ...)
}
\arguments{
  \item{S}{A matrix of the scores of the principal components.}
  \item{k}{A vector containing the sample set sizes to be evaluated.}
  \item{method}{The sampling algorithm. Options are: (i) "kss" (Kennard-Stone Sampling), (ii) "kms" (K-means Sampling) the default, (iii) "clhs" (conditioned Latin hypercube Sampling).}
  \item{repetitions}{The number of times that the sampling must be carried out for each sample size to be evaluated. The results of the final msd is the average of the ones obtained at each iteration. Note that since the "kss" method is deterministic and always return the same results, there is no need for repetitions.}
  \item{n}{The number of equally spaced points at which the probability densities are to be estimated (see \code{\link[stats]{density}}).}
  \item{from, to}{A vector of the left and right-most points of the grid at which the densities are to be estimated. Default is the minimums and maximums  of the variables in `S`.}
  \item{bw}{A vector containing the the smoothing bandwidth to be use for the probability densities (see \code{\link[stats]{density}}).}
  \item{...}{arguments to be passed to the calibration sampling algorithms, i.e. additional aruments to be used for the clhs, kenStone or naes functions which run inside this function}
}
\description{
This function works by comparing the probability density function (pdf) of the population and the pdf of the sample set in order to asses the representativeness of the sample set. See Ramirez-Lopez,  et al. (2014) for more details.}
\references{
Ramirez-Lopez, L., Schmidt, K., Behrens, T., van Wesemael, B., Dematte, J. A., Scholten, T. (2014). Sampling optimal calibration sets in soil infrared spectroscopy. Geoderma, 226, 140-150.
}
\value{
A table with the following columns:
 - css: the sample set size (k)
 - msd
 - msd_sd: the standard deviation of the msd for all the repetitions (does not apply to "kss" since it  always return the same results)

}
