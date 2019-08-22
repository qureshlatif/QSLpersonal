tssr <- function (lat, long, dateTime) {
  require(suncalc)
  require(lutz)
  require(sf)
  # lat and long as decimal degrees with negative value as southern and western
  # date as a character e.g., '2009-04-29'
  # time as a character e.g., '2009-04-29 06:30:00'

  tz <- tz_lookup_coords(lat = lat, lon = long, method = "accurate", 
                         warn = FALSE)
  date <- substr(dateTime, 1, 10)
  tz.unq <- unique(tz)
  tssr <- rep(NA, length(dateTime))
  for(z in tz.unq) {
    ind <- which(tz == z & !is.na(date))
    dat <- data.frame(lat = lat[ind],
                      lon = long[ind],
                      date = as.Date(date[ind], tz = z))
    srt <- getSunlightTimes(data = dat, keep = "sunrise", tz = z)$sunrise
    tssr[ind] <- as.numeric(difftime(as.POSIXlt(dateTime[ind], tz = z), 
                                     srt, units = "mins"))
  }
  return(tssr = tssr)
}