% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateerrorparams.R
\name{estimateErrorParams}
\alias{estimateErrorParams}
\title{Estimate prediction error distribution parameters}
\usage{
estimateErrorParams(
  train_nodes,
  test_nodes,
  mspewhat,
  biaswhat,
  intervalwhat,
  pwhat,
  qwhat,
  mcrwhat,
  alpha,
  n.test
)
}
\arguments{
\item{train_nodes}{A \code{data.table} indicating which out-of-bag prediction
errors are in each terminal node of each tree in the random forest. It
must be formatted like the output of the \code{findOOBErrors} function.}

\item{test_nodes}{A \code{data.table} indicating which test observations are
in each terminal node of each tree in the random forest. It must be
formatted like the output of the \code{findTestPreds} function.}

\item{mspewhat}{A boolean indicating whether to estimate conditional MSPE.}

\item{biaswhat}{A boolean indicating whether to estimate conditional bias.}

\item{intervalwhat}{A boolean indicating whether to estimate conditional
prediction intervals.}

\item{pwhat}{A boolean indicating whether to estimate the conditional
prediction error CDFs.}

\item{qwhat}{A boolean indicating whether to estimate the conditional
prediction error quantile functions.}

\item{mcrwhat}{A boolean indicating whether to estimate the conditional
misclassification rate.}

\item{alpha}{A vector of type-I error rates desired for the conditional prediction
intervals; required if \code{intervalwhat} is \code{TRUE}.}

\item{n.test}{The number of test observations.}
}
\value{
A \code{data.frame} with one or more of the following columns:

  \item{pred}{The random forest predictions of the test observations}
  \item{mspe}{The estimated conditional mean squared prediction errors of
  the random forest predictions}
  \item{bias}{The estimated conditional biases of the random forest
  predictions}
  \item{lower_alpha}{The estimated lower bounds of the conditional alpha-level
  prediction intervals for the test observations}
  \item{upper_alpha}{The estimated upper bounds of the conditional alpha-level
  prediction intervals for the test observations}
  \item{mcr}{The estimated conditional misclassification rate of the random
  forest predictions}

  In addition, one or both of the following functions:

  \item{perror}{The estimated cumulative distribution functions of the
  conditional error distributions associated with the test predictions}
  \item{qerror}{The estimated quantile functions of the conditional error
  distributions associated with the test predictions}
}
\description{
Estimates the prediction error distribution parameters requested in the input
to \code{quantForestError}.
}
\details{
This function is for internal use.
}
\seealso{
\code{\link{quantForestError}}, \code{\link{findOOBErrors}}, \code{\link{findTestPreds}}
}
\author{
Benjamin Lu \code{<b.lu@berkeley.edu>}; Johanna Hardin \code{<jo.hardin@pomona.edu>}
}
\keyword{internal}
