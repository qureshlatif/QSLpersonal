\name{sumJAGS}
\alias{sumJAGS}
\title{
  Provide posterior estimates in list format and basic summaries.
}
\description{
  Takes raw output from saveJAGS package, and provides estimates in list format and a summary table containing means, SDs, quantiles, Rhat, and neffective for each parameter.
}
\usage{
  sumJAGS(rawJAGS, npar.large = 50)
}
\arguments{
  \item{rawJAGS}{Raw output containing posterior estimates in array format produced by saveJAGS::combineSaves.}
  \item{npar.large}{If number of parameters is greater than this (default = 50), Rhat is calculated in chunks to avoid "Error in chol.default(W) :". Additionally, compilation speed for Rhat appears to be optimized at 0 < npar.large < number of parameters, so the user can try setting npar.large at different levels to optimize speed.}
}
\value{A list object with 3 elements: (1) "sims.array" contains raw output as provided by `combineSaves`, (2) "sims.list" contains estimates in list format, and (3) "summary" contains summary table with columns mean, SD, quantiles (0.025, 0.25, 0.5, 0.75, and 0.975), Rhat, and neffective.}
\author{
  Quresh S. Latif, Bird Conservancy of the Rockies
}

\keyword{saveJAGS}
\keyword{JAGS}
