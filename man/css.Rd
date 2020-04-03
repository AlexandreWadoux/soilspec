\name{css}
\alias{css}
\title{Optimal sample size for infrared spectroscopic data}
\usage{
css(S, k, method = "kms", repetitions = 10, n = 512, from, to, bw, ...)
}
\arguments{
  \item{S}{A matrix of the scores of the principal components.}
  \item{k}{A vector containing the sample set sizes to be evaluated.}
  \item{method}{The sampling algorithm. Options are: (i) \code{kss} (Kennard-Stone sampling), (ii) \code{kms} (\emph{K}-means sampling) the default, (iii) \code{clhs} (conditioned Latin Hypercube sampling).}
  \item{repetitions}{The number of times that the sampling must be carried out for each sample size to be evaluated. The results of the final msd is the average of the ones obtained at each iteration. Note that since the \code{kss} method is deterministic and always return the same results, there is no need for repetitions.}
  \item{n}{The number of equally spaced points at which the probability densities are to be estimated (see \code{\link[stats]{density}}).}
  \item{from, to}{A vector of the left and right-most points of the grid at which the densities are to be estimated. Default is the minimums and maximums  of the variables in \code{S}.}
  \item{bw}{A vector containing the the smoothing bandwidth to be use for the probability densities (see \code{\link[stats]{density}}).}
  \item{...}{Arguments to be passed to the calibration sampling algorithms, i.e. additional aruments to be used for the \code{\link[clhs]{clhs}}, \code{\link[prospectr]{kenStone}} or \code{\link[prospectr]{naes}} functions which run inside this function.}
}
\description{
Function to determine the optimal number of spectra to be sent to the laboratory for soil analysis. This function works by comparing the probability density function (pdf) of the population and the pdf of the sample set in order to asses the representativeness of the sample. The two pdfs are compared based on the mean euclidean distance (msd). See Ramirez-Lopez,  et al. (2014) for more details.}
\references{
Ramirez-Lopez, L., Schmidt, K., Behrens, T., van Wesemael, B., Dematte, J. A., Scholten, T. (2014). Sampling optimal calibration sets in soil infrared spectroscopy. Geoderma, 226, 140-150.
}
\value{
A table with the following columns:
  \item{css}{The sample size (\eqn{k}).}
  \item{msd}{Value of the msd for each sample size.}
  \item{msd_sd}{The standard deviation of the msd for all the repetitions (does not apply to \code{kss} since it always return the same results).}

}
