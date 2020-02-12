\name{eval}
\alias{eval}
\title{Series of accuracy measures for quantitative or categorical variables}
\usage{
eval(obs, pred, obj = 'quant', ...)
}
\arguments{
  \item{obs}{A vector of measured values or classes.}
  \item{pred}{A vector of predicted values or classes.}
  \item{obj}{Either `quant` or `cat` to specify whether `obs` and `pred` are categorical or quantitative.}
  \item{...}{additional arguments}
}
\description{
This function returns a set of accuracy measures between observed and predicted values or classes.}
