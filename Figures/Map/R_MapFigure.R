
require(dplyr)
require(sf)
require(ggplot2)
require(ggspatial)
require(viridis)
require(deltamapr)
require(ggrepel)
require(cowplot)

setwd("C:/Users/bmahardja/Documents/GitHub/DeltaSmelt_SummerFallX2_VOI/Figures/Map")

#Add crs for lat and long
crsLONGLAT <- "+proj=longlat +datum=WGS84 +no_defs"

#Read CA state file
CA_shape <- st_read("CA_State.shp") 

#Create the California inset map
fig_inset <- ggplot() + theme_bw()+ geom_sf(data = CA_shape, fill="white", color="black") +
  coord_sf(crs=crsLONGLAT)+
  geom_rect(aes(xmin = -122.5, xmax = -121.4, ymin = 37.65, ymax = 38.61), fill = NA, color = "blue", size = 1) +
  annotate(geom = "text", x = -118.898056, y = 36.154072, label="California",size=8,angle = 315) + 
  theme_void() +
  theme(panel.background = element_rect(fill = "gray80", color = NA))
fig_inset



#Grab X2 data from deltmapr package
P_X2_data = subset(P_X2, RKI/5 == round(RKI/5)) %>% filter(RKI<=100)

# Create a linestring from the points
line_coords <- st_coordinates(P_X2_data)
line_sf <- st_sf(geometry = st_sfc(st_linestring(line_coords), crs = 4326))

#Read ocean boundary shape file to extend Golden Gate Bridge area
Ocean<-st_read("GIS_ADMIN_OCEAN_BAY.shp") %>% filter(LANDNAME=="Pacific Ocean")

#Create polygon for Suisun Bay and Marsh
SuisunMarsh <- R_EDSM_Subregions_Mahardja_FLOAT %>% filter(SubRegion %in% c("Suisun Marsh"))
test<-R_EDSM_Subregions_Mahardja_FLOAT
SuisunBay <- R_EDSM_Subregions_Mahardja_FLOAT %>% filter(SubRegion %in% c("Mid Suisun Bay", "Grizzly Bay","Honker Bay","West Suisun Bay")) %>% st_union()
  
Delta <- R_Delta
#Create label for certain regions
DT = data.frame(
  lat=c(38.201114, 38.066950, 38.046551,37.81934,38.187636),
  long=c(-121.993206, -122.007843, -121.604598,-122.4778,-121.975941),
  name=c("Suisun Marsh","Suisun Bay","Sacramento-San Joaquin Delta","Golden Gate Bridge","Belden's Landing")
)

DT = st_as_sf(DT, coords = c("long","lat"), remove = FALSE,crs=crsLONGLAT)

PT<-data.frame(
  lat=c(38.187636),
  long=c(-121.975941),
  name=c("Belden's Landing")
)
PT = st_as_sf(PT, coords = c("long","lat"), remove = FALSE,crs=crsLONGLAT)


fig1<-ggplot() + theme_bw()+
  geom_sf(data = Ocean, fill="deepskyblue", color="deepskyblue") +
  geom_sf(data = SuisunMarsh, alpha=0.5,color="black",fill="black",size=2) +
  geom_sf(data = Delta, alpha=0.5,color="grey",fill="grey",size=2) +
  geom_sf(data = SuisunBay, alpha=0.5,color="blue",fill="blue",size=2) +
  geom_sf(data = WW_Delta, fill="deepskyblue", color="deepskyblue") +
  #geom_sf(data = R_DSIBM, fill=NA, color="grey",alpha=0.1, size=2) +
  geom_sf(data=P_X2_data)+
  geom_sf(data=line_sf)+
  geom_sf(data=PT)+
  geom_sf_text(data = P_X2_data, 
               aes( label = RKI),vjust = 1.4,hjust=0.3)+
  coord_sf(xlim = c(-122.5, -121.4), ylim = c(37.65, 38.61),crs=crsLONGLAT)+
  geom_label_repel(data=DT, aes(x=long,y=lat,label=name), nudge_x=c(0.1,0.1,0.02,0.2), nudge_y=c(0.1,-0.1,-0.2,-0.12)
                   ,segment.alpha=0.7,color="black", size=3.7,segment.linetype="dashed")  +
  annotation_north_arrow(location = "tr", which_north = "true", pad_y = unit(1.0, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "tr", width_hint = 0.5)+  
  theme(axis.text.x = element_text(size=12, color="black"),axis.text.y = element_text(size=12, color="black"),
       axis.title.x=element_blank(),axis.title.y=element_blank())
fig1


# Combine main and inset maps
combined_map <- ggdraw() +
  draw_plot(fig1) +
  draw_plot(fig_inset, x = 0.15, y = 0.65, width = 0.3, height = 0.3)  # Adjust position and size of inset
combined_map


#Print out the map
tiff(filename="Figure01_Map.tiff", units="in",type="cairo", bg="white", height=10, 
     width=11, res=300, compression="lzw")
combined_map
dev.off()



