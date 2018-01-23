# 5-1-2016 MRC-Epid JHZ

l <- read.csv(l,file="map.csv")

library(maps)
library(mapdata)
data(world.cites)
country <- c("France","Italy","Spain","UK","Netherlands","Greece","Germany","Sweden","Denmark")
caps <- subset(world.cities,country.etc%in%country&capital==1)

with(l[-(25:26),],{
  png("map.png",width=7,height=7,res=300,units="in")
  map('world', country)
# map.text('world',country)
  with(caps, {
    points(long,lat,col=14,pch=20)
    text(long,lat+0.5,name,cex=0.5,col="red")
  })
  text(lon,lat,letters[-(25:26)],cex=0.5,col="blue")
  xlim <- range(lon)
  ylim <- range(lat)
  dev.off()
})
with(l[-(25:26),],paste(letters[-(25:26)],"-",centre))

