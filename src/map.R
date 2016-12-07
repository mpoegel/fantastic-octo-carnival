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

pal <- colorRampPalette(c("light gray", "red"))


#
# Plot Minnesota maps
#

mn.ground.truth <- read.csv("data/external/MN-cities.csv", header=FALSE,
                            col.names = c("lat", "lon", "name"))

mn.map.data <- map_data("state", region = "minnesota")

p1 <- ggplot() +
  geom_map(data = mn.map.data, map = state.map.data,
           aes(x=long, y=lat, map_id=region), fill = light.gray, color = dark.gray) +
  geom_point(data = mn.ground.truth, stat = "identity", aes(x = lon, y = lat), color = "blue") +
  geom_label_repel(data = mn.ground.truth, aes(x = lon, y = lat, label = name), size=3,
                   box.padding = unit(0.1, "cm")) +
  ggtitle("Minnesota - Sample Cities") +
  mapTheme()
p1

# Plot the ground truth against the output
mn.ground.truth <- read.table("data/processed/Minnesota-latlong.dat", header=FALSE, skip=1, sep=" ",
                            col.names = c("id", "lat", "lon"))
mn.estimates <- read.csv("output/minnesota.yhat")

# png("output/mn_ground_truth_vs_estimates.png")
p2 <- ggplot() +
  geom_map(data = mn.map.data, map = state.map.data,
           aes(x=long, y=lat, map_id=region), fill = light.gray, color = dark.gray) +
  geom_text(data = mn.ground.truth, aes(x = lon, y = lat, label = id), size = 3, color = "blue") +
  geom_text(data = mn.estimates, aes(x = lon, y = lat, label = id), size = 3, color = "red") +
  ggtitle("Minnesota - Ground Truth v. Estimate") +
  mapTheme()
p2
# dev.off()

# Plot the centrality index
mn.ci <- read.csv("output/minnesota.ci")
colors <- pal(length(mn.ci))

p3 <- ggplot() +
  geom_map(data = mi.map.data, map = state.map.data,
           aes(x=long, y=lat, map_id=region), fill = light.gray, color = dark.gray) +
  geom_point(data = mn.ci, stat = "identity", aes(x = lon, y = lat, color = CI)) +
  scale_colour_gradient(limits=c(min(mn.ci$CI), max(mn.ci$CI)), high="red") +
  ggtitle("Minnesota - Centrality Index") +
  mapTheme()
p3

# COARSENED: Plot the ground truth against the output
mn.estimates <- read.csv("output/minnesota.coarse.yhat")

# png("output/mn_ground_truth_vs_estimates.png")
p4 <- ggplot() +
  geom_map(data = mn.map.data, map = state.map.data,
           aes(x=long, y=lat, map_id=region), fill = light.gray, color = dark.gray) +
  geom_text(data = mn.ground.truth, aes(x = lon, y = lat, label = id), size = 3, color = "blue") +
  geom_text(data = mn.estimates, aes(x = lon, y = lat, label = id), size = 3, color = "red") +
  ggtitle("Minnesota - Ground Truth v. Estimate (Coarsened)") +
  mapTheme()
p4
# dev.off()

# COARSENED: Plot the centrality index
mn.ci <- read.csv("output/minnesota.coarse.ci")
colors <- pal(length(mn.ci))

p5 <- ggplot() +
  geom_map(data = mi.map.data, map = state.map.data,
           aes(x=long, y=lat, map_id=region), fill = light.gray, color = dark.gray) +
  geom_point(data = mn.ci, stat = "identity", aes(x = lon, y = lat, color = CI)) +
  scale_colour_gradient(limits=c(min(mn.ci$CI), max(mn.ci$CI)), high="red") +
  ggtitle("Minnesota - Centrality Index (Coarsened)") +
  mapTheme()
p5
