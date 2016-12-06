library(extrafont)
library(ggplot2)
library(ggrepel)
library(grid)
library(maps)
library(mapdata)
library(ggmap)


mapTheme <- function() {
  theme(
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )
}


#
# Plot Minnesota maps
#

mi.ground.truth <- read.csv("data/external/MN-cities.csv", header=FALSE,
                            col.names = c("lat", "lon", "name"))

mi.map.data <- map_data("state", region = "minnesota")

p1 <- ggplot() +
  geom_map(data = mi.map.data, map = state.map.data,
           aes(x=long, y=lat, map_id=region), fill = light.gray, color = dark.gray) +
  geom_point(data = mi.ground.truth, stat = "identity", aes(x = lon, y = lat), color = "blue") +
  geom_label_repel(data = mi.ground.truth, aes(x = lon, y = lat, label = name), size=3,
                   box.padding = unit(0.1, "cm")) +
  mapTheme()
p1

# Plot the ground truth against the output
mi.ground.truth <- read.table("data/processed/Minnesota-latlong.dat", header=FALSE, skip=1, sep=" ",
                            col.names = c("id", "lat", "lon"))
mi.estimates <- read.csv("output/mn-yhat.dat")

p2 <- ggplot() +
  geom_map(data = mi.map.data, map = state.map.data,
           aes(x=long, y=lat, map_id=region), fill = light.gray, color = dark.gray) +
  geom_text(data = mi.ground.truth, aes(x = lon, y = lat, label = id), size = 3, color = "blue") +
  geom_text(data = mi.estimates, aes(x = lon, y = lat, label = id), size = 3, color = "red") +
  ggtitle("Minnesota - Ground Truth v. Estimate") +
  mapTheme()
p2

