# convert units
`%#%` <- function(x,y)
{
  # convert to si units
  if(is.numeric(x))
  {
    num <- x
    name <- y
    pow <- +1
  }
  else # convert from si units
  {
    if(!is.numeric(y)) { return(x %#% 1 %#% y) }
    
    num <- y
    name <- x
    pow <- -1
  }
  
  name <- ctmm:::canonical.name(name)
  
  # DIV <- grepl("/",name)
  # if(DIV)
  # {
  #   #
  # }
  
  OP <- getOption("time.units")
  
  alias <- list()
  scale <- c()
  
  add <- function(a,s)
  {
    n <- length(alias)
    alias[[n+1]] <<- ctmm:::canonical.name(a)
    scale[n+1] <<- s
  }
  
  # this database should be generated on load instead of on function call
  
  # TIME
  add(c("\u03BCs","\u03BCs.","microsecond","microseconds"),1E-6)
  add(c("ms","ms.","milisecond","miliseconds"),1/1000)
  add(c("s","s.","sec","sec.","second","seconds"),1)
  add(c("min","min.","minute","minutes"),60)
  add(c("h","h.","hr","hr.","hour","hours"),60^2)
  add(c("day","days"),UNITS[[OP]]$DAY) #day
  add(c("wk","wk.","week","weeks"),7*UNITS[[OP]]$DAY) # week
  add(c("mon","mon.","month","months"),UNITS[[OP]]$MONTH) # month
  add(c("yr","yr.","year","years"),UNITS[[OP]]$YEAR) # year
  add(c("ka","ky","millennium","millenniums","millennia","kiloannum","kiloannums","kiloyear","kiloyears"),1000*UNITS[[OP]]$YEAR)
  add(c("ma","my","megaannum","megaannums","megaanna","megayear","megayears","millionennium","millionenniums","millionennia"),1000^2*UNITS[[OP]]$YEAR)
  add(c("ae","ga","gy","gyr","aeon","aeons","eon","eons","gigayear","gigayears","gigaannum","gigaannums","giggaanna"),1000^3*UNITS[[OP]]$YEAR)
  
  # Distance conversions
  add(c("\u03BCm","\u03BCm.","micron","microns","micrometer","micrometers"),1E-6)
  add(c("mm","mm.","milimeter","milimeters"),1/1000)
  add(c("cm","cm.","centimeter","centimeters"),1/100)
  add(c("m","m.","meter","meters"),1)
  add(c("km","km.","kilometer","kilometers"),1000)
  add(c("in","in.","inch","inches"),0.3048/12)
  add(c("ft","ft.","foot","feet"),0.3048)
  add(c("yd","yd.","yard","yards"),0.3048*3)
  add(c("mi","mi.","mile","miles"),0.3048*5280)
  
  # Area conversions
  add(c("\u03BCm\u00B2","\u03BCm.\u00B2","micron\u00B2","microns\u00B2","micrometer\u00B2","micrometers\u00B2","\u03BCm^2","\u03BCm.^2","micron^2","microns^2","micrometer^2","micrometers^2","square micron","square microns","square micrometer","square micrometers","micron squared","microns squared","micrometer squared","micrometers squared"),1E-12)
  add(c("mm\u00B2","mm.\u00B2","milimeter\u00B2","milimeters\u00B2","mm^2","mm.^2","milimeter^2","milimeters^2","square milimeter","square milimeters","milimeter squared","milimeters squared"),1/1000^2)
  add(c("cm\u00B2","cm.\u00B2","centimeter\u00B2","centimeters\u00B2","cm^2","cm.^2","centimeter^2","centimeters^2","square centimeter","square centimeters","centimeter squared","centimeters squared"),1/100^2)
  add(c("m\u00B2","m.\u00B2","meter\u00B2","meters\u00B2","m^2","m.^2","meter^2","meters^2","square meter","square meters","meter squared","meters squared"),1)
  add(c("ha","hectare","hectares","hm\u00B2","hectometer\u00B2","hectometre\u00B2","hectometers\u00B2","hectometres\u00B2","hm^2","hectometer^2","hectometre^2","hectometers^2","hectometres^2","square hm","square hectometer","square hectometre","square hectometers","square hectometres"),100^2)
  add(c("km\u00B2","km.\u00B2","kilometer\u00B2","kilometers\u00B2","km^2","km.^2","kilometer^2","kilometers^2","square kilometer","square kilometers","kilometer squared","kilometers squared"),1000^2)
  add(c("in\u00B2","in.\u00B2","inch\u00B2","inches\u00B2","in^2","in.^2","inch^2","inches^2"),(0.3048/12)^2)
  add(c("ft\u00B2","ft.\u00B2","foot\u00B2","feet\u00B2","ft^2","ft.^2","foot^2","feet^2","square foot","square feet","foot squared","feet squared"),0.3048^2)
  add(c("yd\u00B2","yd.\u00B2","yard\u00B2","yards\u00B2","yd^2","yd.^2","yard^2","yards^2","square yard","square yards","yard squared","yards squared"),(0.3048*3)^2)
  add(c("mi\u00B2","mi.\u00B2","mile\u00B2","miles\u00B2","mi^2","mi.^2","mile^2","miles^2","square mile","square miles","mile squared","miles squared"),(0.3048*5280)^2)
  
  # speed
  add(c("mps","m/s","m/sec","meter/sec","meter/second","meters/second"),1)
  add(c("kmph","kph","km/h","km/h","km/hr","kilometer/hour","kilometers/hour"),0.277777777777777777777)
  add(c("mph","mi/h","mi/hr","mile/h","mile/hr","mile/hour","miles/hour"),0.44704)
  add(c("fps","ft/s","ft/sec","feet/second"),0.3048)
  add(c('kt','kn','knot','knots'),1.852 * 0.277777777777777777777)
  add(c("km/s","kmps","km/sec"),1000)
  add(c("cm/s","cmps","cm/sec"),1/100)
  add(c("mm/s","mmps","mm/sec"),1/1000)
  add(c("\u03BCm/s","\u03BCmps","\u03BCm/sec"),1E-6)
  
  for(i in 1:length(alias))
  {
    if(name %in% alias[[i]]) { return(num*scale[i]^pow) }
  }
  stop(paste("Unit",name,"unknown."))
}