library(ape)
library(geiger)
library(dplyr)
library(sf)
library(ggplot2)
library(picante)
library(rnaturalearth)
library(rgdal)
library(biscale)
library(cowplot)


Tree <- read.nexus("4000phylogeny.tre") ##load in tree

Occ <- read.csv("madagascar_ferns.csv") ##load in occurrences

Occ <- st_as_sf(Occ, coords=c("decimalLongitude", "decimalLatitude"), ##convert occurrences to sf object
                crs='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')

behrmann <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" ##set projection for reprojection

Occ <- st_transform(Occ, behrmann) ##reproject occurrence data

Map <- ne_countries(scale = 50, returnclass = "sf", ##obtain map
                    type = "countries") %>% 
  filter(sovereignt == "Madagascar")

Map <- st_transform(Map, behrmann) ##reproject map

occMap <-   ##plot occurrence map
  ggplot() +
  geom_sf(data = Map, color = "black", size = 0.05, fill = "#f7f7f7")+
  geom_sf(data = Occ , aes(geometry = geometry), size = 0.25) +
  geom_sf(data = Map, color = "black", size = 0.05, fill = NA)+
  labs(x="Longitude", y="Latitude", color = "Species Group")+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme(panel.grid.major = element_line(colour = "#c9c9c9", 
                                        linetype = "dashed", 
                                        size = 0.2), 
        panel.background = element_rect(fill = "#f0f8ff"), 
        panel.border = element_rect(fill = NA),
        axis.text.x = element_text(angle = 305, vjust = 0.2, hjust=0.2),
        legend.key = element_blank())
occMap


MatchedNames <- intersect(Tree$tip.label, Occ$scientificName) ##match names in tree and occurrences

Occ <- Occ %>% filter(scientificName %in% MatchedNames) ##retain occurrences of only species in tree

GridSize <- 50000 ##specify grid dimensions (in this case, in meters)

Grid <- st_make_grid(x = Map, what = "polygons", cellsize = GridSize) ##make grid

Grid <- st_sf(idcell = 1:length(Grid), geom = Grid) %>%  ##convert grid to sf object
  st_cast("POLYGON")

IntersectedGrid <- st_intersection(Grid,Occ)  ##intersect grid and occurrences

IntersectedGridDF <- as.data.frame(IntersectedGrid)[,c(1,2)] 

GridSpeciesMatrix <- merge.data.frame(x = data.frame(idcell = Grid$idcell), ##create species community matrix
                                      y = IntersectedGridDF, 
                                      all.x = TRUE, by.x = TRUE, 
                                      sort = TRUE) %>% table()

PdOutput <- pd(GridSpeciesMatrix,Tree) ##calculate PD

PdOutput <- PdOutput %>%  ##assign cell IDs
  mutate(idcell = row_number())

PdOutput$idcell <- as.integer(PdOutput$idcell)

GridPD <- left_join(Grid,PdOutput, by = "idcell") ##join PD data with grid

ClippedGridPD <- st_intersection(Map, GridPD) ##drop cells outside of study area

ClippedGridPD <- ClippedGridPD %>% select(idcell, PD, SR, geometry)

lims_map <- st_buffer(Map, dist = 100) %>% st_bbox() ##define extent of map for plotting

gridPD <-   ##plot phylogenetic diversity map
  ggplot() +
  geom_sf(data = Map, color = "black", size = 0.05)+
  geom_sf(data=ClippedGridPD, aes(fill = PD), color = NA) +
  geom_sf(data = Map, color = "black", size = 0.2, fill=NA)+
  scale_fill_distiller(palette="Reds", direction = 1, limit=c(1,max(GridPD$PD)),values=c(0,1), na.value = "#ffffff")+
  coord_sf(
    xlim = c(lims_map["xmin"], lims_map["xmax"]),
    ylim = c(lims_map["ymin"], lims_map["ymax"])
  ) +
  scale_x_continuous(breaks = c(-84)) +
  theme(
    plot.background = element_rect(fill = "#ffffff"),
    panel.background = element_rect(fill = "#f0f8ff"),
    panel.grid = element_blank(),
    line = element_blank(),
    rect = element_blank()
  ) + labs(fill = "PD (Ma)")
gridPD

gridSR <- ##plot species richness map
  ggplot() +
  geom_sf(data = Map, color = "black", size = 0.05)+
  geom_sf(data=ClippedGridPD, aes(fill = SR), color = NA) +
  geom_sf(data = Map, color = "black", size = 0.2, fill=NA)+
  #scale_fill_scico(direction = 1, palette="vikO", begin=0.5, limits=c(1,max(madGridPD$SR)),na.value = "#ffffff") +
  scale_fill_distiller(palette="Reds", direction = 1, limit=c(1,max(GridPD$SR)),values=c(0,1), na.value = "#ffffff")+
  coord_sf(
    xlim = c(lims_map["xmin"], lims_map["xmax"]),
    ylim = c(lims_map["ymin"], lims_map["ymax"])
  ) +
  scale_x_continuous(breaks = c(-84)) +
  theme(
    plot.background = element_rect(fill = "#ffffff"),
    panel.background = element_rect(fill = "#f0f8ff"),
    panel.grid = element_blank(),
    line = element_blank(),
    rect = element_blank()
  ) + labs(fill = "Species Richness")
gridSR


BiscaleData <- bi_class(ClippedGridPD, x = PD, y = SR,   ##break data into quantiles
                        style = "fisher", dim = 3)


Legend <- bi_legend(pal = "GrPink",     ##make bivariate legend for plotting
                    dim = 3,
                    xlab = "Higher PD ",
                    ylab = "Higher SR",
                    size = 8)

BivariateMap <- ggplot() +      ##plot map
  geom_sf(data = BiscaleData, mapping = aes(fill = bi_class), color = "white", size = 0.1, show.legend = F) +
  bi_scale_fill(pal = "GrPink", dim = 3)+
  bi_theme()


FinalPlot <- ggdraw() +     ##finalize map & legend
  draw_plot(BivariateMap, 0, 0, 1, 1) +
  draw_plot(Legend, 0.55, 0.1, 0.275, 0.275)
FinalPlot