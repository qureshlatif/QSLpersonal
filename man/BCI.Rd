\name{BCI}
\alias{BCI}
\title{Summarize Bayesian posterior parameter distribution.}
\description{Summarizes posterior samples as median and credible intervals using specified percentiles.}
\usage{
BCI(x, ndig = 2, BCIpercent = 95, flag.sig = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Numeric vector of samples drawn from posterior parameter distribution.}
  \item{ndig}{Number of significant digits to display in summary.}
  \item{BCIpercent}{Percent of posterior samples to be included in credible interval.}
  \item{flag.sig}{If 'TRUE' asterisk is added to summaries with BCIs that exclude zero.}
}
\value{Character string: "[median] ([lower bound], [upper bound])".}
\author{
Quresh S. Latif
}
