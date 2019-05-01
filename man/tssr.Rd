\name{tssr}
\alias{tssr}
\title{Calculate time since sunrise.}
\description{Looks up sunrise for a given latitude, longditude, and date, and then converts 24-hr time to time since sunrise in minutes.}
\usage{
tssr(lat,long,time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lat}{latitude as decimal degrees with negative values representing southern locations.}
  \item{long}{longitude as decimal degrees with negative values representing western locations.}
  \item{dateTime}{date and time as a character e.g., '2009-04-29 06:30:00'.}
}
\value{A value or vector matching the dimensions of dateTime representing the number of minutes from sunrise.}
\author{
Adam Green, Quresh Latif
}
