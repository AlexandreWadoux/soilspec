\name{eval}
\alias{eval}
\title{Accuracy measures for quantitative or categorical variables}
\usage{
eval(obs, pred, obj, ...)
}
\arguments{
  \item{obs}{A vector of measured values or classes.}
  \item{pred}{A vector of predicted values or classes.}
  \item{obj}{Either \code{quant} or \code{cat} to specify whether \code{obs} and \code{pred} are quantitative or categorical.}
  \item{...}{Further not yet used parameters.}
}
\description{
Computation of a set of accuracy measures between observed and predicted values or classes.
}
\details{

The accuracy measures for quantitative variables (when \code{obj = 'quant'}) are
\itemize{
\item{ME: mean error}
\item{RMSE: root mean square error}
\item{r2: squared Pearson correlation coefficient}
\item{R2: coefficient of determination}
\item{rhoC: concordance correlation coefficient}
\item{RPD: ratio of performance to deviation}
\item{RPIQ: ratio of performance to inter-quartile distance}
}
The accuracy measures for categorical variables (when \code{obj = 'cat'}):
\itemize{
\item{OA: overall accuracy}
\item{kappa: Cohen's kappa statistic}
}
}

\value{
Data frame with accuracy measures.
}
