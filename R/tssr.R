tssr <- function(lat,long,dateTime) {
  require(suncalc)
  require(lutz)
  require(sf)
  
  # lat and long as decimal degrees with negative value as southern and western
  # date as a character e.g., '2009-04-29'
  # time as a character e.g., '2009-04-29 06:30:00'
  tz <- tz_lookup_coords(lat=lat,lon=long,method='accurate',warn=FALSE)
  date <- substr(dateTime, 1, 10)
  srt <- getSunlightTimes(date=as.Date(date,tz=tz),lat=lat,lon=long,keep='sunrise',tz=tz)$sunrise
  
  tssr <- as.numeric(difftime(as.POSIXlt(dateTime,tz=tz), srt, units = 'mins'))
  return(tssr=tssr)
}
