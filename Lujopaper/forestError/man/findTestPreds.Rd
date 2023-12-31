% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findtestpreds.R
\name{findTestPreds}
\alias{findTestPreds}
\title{Compute and locate test predictions}
\usage{
findTestPreds(forest, X.test, n.cores = 1)
}
\arguments{
\item{forest}{The random forest object being used for prediction.}

\item{X.test}{A \code{matrix} or \code{data.frame} with the observations to
be predicted. Each row should be an observation, and each column should be
a predictor variable.}

\item{n.cores}{Number of cores to use (for parallel computation in \code{ranger}).}
}
\value{
A \code{data.table} with the following four columns:

  \item{rowid_test}{The row ID of the test observation as provided by \code{X.test}}
  \item{pred}{The random forest prediction of the test observation}
  \item{tree}{The ID of the tree of the random forest}
  \item{terminal_node}{The ID of the terminal node of the tree in which the
  test observation falls}
}
\description{
Predicts each test observation's response using the random forest and, for
each test observation and tree, finds the terminal node of the tree in which
the test observation falls.
}
\details{
This function accepts regression random forests built using the \code{randomForest},
\code{ranger}, \code{randomForestSRC}, and \code{quantregForest} packages.
}
\seealso{
\code{\link{findOOBErrors}}, \code{\link{quantForestError}}
}
\author{
Benjamin Lu \code{<b.lu@berkeley.edu>}; Johanna Hardin \code{<jo.hardin@pomona.edu>}
}
\keyword{internal}
