\name{VIF}
\alias{VIF}
\title{Calculate variance inflation factors for assessing multicollinearity.}
\description{Calculates variance inflation factors for all indicated variables appearing in data frame.}
\usage{
VIF(data, vars = names(data), returnModels = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data frame containing covariates of interest for multiple regression or multivariate regression-based analysis.}
  \item{vars}{Character vector listing names of covariates in 'data' for which variance inflation factors are desired. Default is all columns in 'data'.}
  \item{returnModels}{If 'TRUE' model objects regressing each variable against all other variables in 'vars' are returned in list form. If 'FALSE' (default), only VIF values are returned.}
}
\value{If returnModels == TRUE, a list whose first element ('vif') is a vector of length = length(vars) containing variance inflation factors (VIF = 1 / (1 - R^2)), and second element is a list of length = length(vars) with lm model objects as elements. If returnModels == FALSE, only VIFs are returned as numeric vector.}
\author{
Quresh S. Latif
}
