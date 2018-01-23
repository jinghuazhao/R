#5-1-2016 MRC-Epid JHZ

library(foreign)
l <- read.csv("map.csv")
library(ggmap)
map <- get_map(location = 'Europe', zoom = 4)
cp <- ggmap(map)+geom_point(aes(x=lon, y=lat, size=N), data=l, alpha=.5, color="red")
png("ggmap.png",height=7,width=7,units="in",res=300)
plot(cp)
dev.off()

