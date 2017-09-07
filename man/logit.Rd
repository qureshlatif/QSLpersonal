\name{logit}
\alias{logit}
\title{Perform logit transformation.}
\description{Applies logit transformation to one or more numbers with range 0-1, resulting in equivalent set of unbounded real numbers with logistic (S-shaped) relationship to inputs.}
\usage{
logit(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Any numeric value, vector, or array whose elements, x, such that 0 < x < 1.}
}
\value{A value, vector, or array matching the dimensions of x and with elements, y, such that -Inf < y < Inf.}
\author{
Quresh S. Latif
}

\seealso{expit}
