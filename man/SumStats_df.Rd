\name{SumStats_df}
\alias{SumStats_df}
\title{Tabulate descriptive statistics.}
\description{Tabulates descriptive statistics for all indicated variables in a data frame. Continuous covariates are summarized as "Mean (SD, range)". Binary covariates are summarized as the percentage of 1s.}
\usage{
SumStats_df(data, vars = names(data), binary = c())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data frame containing covariates of interest for compiling descriptive statistics.}
  \item{vars}{Character vector listing names of covariates in 'data' for which descriptive statistics are desired. Default is all columns in 'data'.}
  \item{binary}{Integer vector indexing elements in 'vars' that are binary covariates. Default is an empty vector.}
}
\value{Matrix with nrow = length(vars) and 2 columns. The first column labeled "Mean (SD, range)" contains summary statistics. The second column contains the sample size for calculating summary statistics (number of non-missing values). Summary statistics for variables indicated in 'binary' will be a percent of 1s.}
\author{
Quresh S. Latif
}
