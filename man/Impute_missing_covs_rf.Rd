\name{Impute_missing_covs_rf}
\alias{Impute_missing_covs_rf}
\title{Impute missing covariate values.}
\description{Imputes missing covariate values using a Random Forest model informed with specified covariates.}
\usage{
Impute_missing_covs_rf(dat, v.fill, v.inform)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{Data frame or matrix containing covariate values.}
  \item{v.fill}{String naming the variable in 'dat' with missing values to be imputed.}
  \item{v.inform}{String vector naming the covariates in 'dat' to be used for informing a random forest model to impute missing values.}
}
\value{Data frame or matrix of same dimensions as 'dat' but with missing values imputed.}
\author{
Quresh S. Latif
}
\note{R packages required: randomForest, dplyr, stringr}
