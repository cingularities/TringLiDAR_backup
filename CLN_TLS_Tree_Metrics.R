# Script segments a point cloud into clusters using watershed segmentation on rasterized point cloud
#remotes::install_github("Jean-Romain/lidR")
#remotes::install_github('tiagodc/TreeLS')
#install.packages(c("dplR","dplyr","tidr", "geosphere", "lidR", "raster", "TreeLS", "rgdal", "rgeos", "sp","sf"))
Packages <- c("dplR","dplyr","tidr", "geosphere", "lidR", "raster", "TreeLS", "rgdal", "rgeos", "sp","sf")
lapply(Packages, library, character.only = TRUE)

library(dplR)
library(sp)
library(EBImage)
library(tidyverse)
library(raster)
library(rgdal) 
library(lidR) 
library(ForestTools)
library(rLiDAR)
library(rTLS)
library(rgl)
library(ggplot2)
library(readxl)
library(geosphere)
library(lidaRtRee)
library(rgdal)
library(tidyverse)
library(sf)
library(raster)
library(lidR)
library(TreeLS)


setwd("D:/projects/Babst_TLidar_Tree_Rings/CLNorton")
packageVersion("TreeLS")
packageVersion("lidR")
packageVersion("raster")
packageVersion("rgdal")
packageVersion("sf")
packageVersion("rgeos")
packageVersion("sp")

library(raster)
library(lidR)
library(TreeLS)
library(rgdal) 
library(sf) 
library(rgeos) 
library(sp)


############################################ Tree Ring and Ground Fusion ####################################################

#select the path for folder where you have the file to load
ground <- read_excel("D:/projects/Babst_TLidar_Tree_Rings/CLNorton/Europe TLS Data/Serial_sampling/metadata_hainich_Mar2019.xlsx")
tree<-read.rwl("D:/projects/Babst_TLidar_Tree_Rings/CLNorton/Europe TLS Data/TLS Data/Hainich_Treerings/Good_FASY.rwl") 
points = lidR::readLAS("D:/projects/Babst_TLidar_Tree_Rings/TLS_SantaRita_UBplot_pointcloud_041922_aligned_CC.las")



###FUNCTION 1 - detrending and chronology###
#check the class, you should see "rwl  "data.frame"
class(tree)
#some basic statistics
dim(tree)
rwl.stats(tree)
options(max.print=1000000)
#Spline
detrend.rwi <- detrend(rwl = tree,method = "Spline")
head(detrend.rwi)
data.crn <- chron(detrend.rwi, prefix = "Chr", prewhiten = TRUE)






###FUNCTION 2 - tree ring data parse and formatting###
#raw tree ring data parsing and mean tree cores
tree.data.parse <- detrend.rwi %>% 
  t() %>%
  as.data.frame() %>%
  rename_with( ~ paste0("year_", .x)) %>%
  rownames_to_column()%>%
  separate(rowname,
           into = c("site", "coreID"),
           sep = "(?<=[A-ZA-Z])(?=[0-9])") %>%
  separate(coreID,
           into = c("treeID", "core"),
           sep = "(?<=[0-9])(?=[A-ZA-Za-z])")


#write.csv(data.transpose, file = "datatranspose.csv")






##RAW GROUND DATA
coor <- c(51.079167, 10.453)
tree23 <- destPoint(coor,110, 76)
plotcenter <- destPoint(tree23,98, 28.1)
MOAB_25 <- c(-109.344632,38.26229)
ground_moab <- read.csv("D:/projects/RaBET/NEON_MOAB/ground_data.csv")
###FUNCTION 3 - ground data azimuth and distance to coordinates###
tls_ground_coor <- function(coor, dist_azim) {
  
  coor_list <- list() #creates empty list for new coordinates
  
  for (i in 1:nrow(dist_azim)) {
    stemAzimuth <- dist_azim$Azimuth[i]
    stemDistance <- dist_azim$Distance[i]
    
    new_coordinates <- destPoint(coor, stemAzimuth,  stemDistance) %>%
      as.data.frame() %>% 
      add_column(treeID = dist_azim$TreeID[i]) %>%
      add_column(stemAzimuth = stemAzimuth) %>%
      add_column(stemDistance = stemDistance) %>%
      add_column(species = dist_azim$Species[i]) %>%
      add_column(Height = dist_azim$Height[i])%>%
      add_column(DBH = dist_azim$DBH[i]) 
    
    
    coor_list[[i]] <- new_coordinates #put result of loop in list
    
  }
  
  coor_df<- coor_list %>% bind_rows() %>% na.omit() 
  
  
  return(coor_df)
}

Hainich_plot <- tls_ground_coor(coor, ground) 
moab_plot <- tls_ground_coor(MOAB_25,ground_moab)
write.csv(UTM_moab, file = "ground_validation_moab_utm.csv")
setwd("D:/projects/RaBET/NEON_MOAB/")

###tree_ring_ground_spatialpoints 1
xy_ground <- data.frame(X = moab_plot$lon, Y = moab_plot$lat) ##CONVERTLATTOLONG
coordinates(xy_ground) <- ~ X + Y
proj4string(xy_ground) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  ## for example
tree.points.ground <- spTransform(xy_ground, CRS("+proj=utm +zone=12 +ellps=WGS84"))
UTM_moab <- tree.points.ground@coords %>% as.data.frame()


###FUNCTION 4 -  tree ring and ground data merge###
###TREE rings and ground merge#####
treering_ground <- tree.data.parse %>%
  subset(site == "HF") %>%
  lapply(as.numeric)%>%
  as.data.frame() %>%
  group_by(treeID) %>%
  summarise(across(-c(site, core), mean, na.rm = TRUE)) %>%
  rev()%>%
  left_join(Hainich_plot, by = "treeID") 









###FUNCTION 5 - DBH and Percentage reconstruction###
###DBH RECONSTRUCITONS
diameter.recon.ground.tree <- function(treering_ground){
  treering_ground[is.na(treering_ground)] <- 0
  ###DBH### 
  # Identify the columns that have data of interest
  year_cols <- which(substr(x = colnames(treering_ground), start = 1, stop = 4) == "year")
  # Add column with total
  treering_ground$lp <- rowSums(x = treering_ground[, year_cols], na.rm = TRUE)
  # Make a copy of that data frame that we will use for percentages
  lh <- treering_ground
  lh[, year_cols] <- sapply(1:ncol(lh[, year_cols]), function(col){rowSums(lh[, year_cols][1:col])})
  diff <- lh
  # Do percentage calculations for only those columns with data of interest
  diff[, year_cols] <- (diff$lp) - (diff[, year_cols])
  fraction <- lh
  fraction[, year_cols] <- diff[, year_cols]/diff$lp
  diameter <- lh
  diameter[, year_cols] <- as.numeric(fraction$DBH)*fraction[, year_cols]

  # Reality check to see that percentages add up to 100
  
return(diameter)
}
  

DBH_estimate <- diameter.recon.ground.tree(treering_ground)


##PERCENTAGE DBH####
DBH_perc_es
# Identify the columns that have data of interest
year_cols <- which(substr(x = colnames(DBH_estimate), start = 1, stop = 4) == "year")
# Make a copy of that data frame that we will use for percentages
data.perc <- DBH_estimate
# Do percentage calculations for only those columns with data of interest
data.perc[, year_cols] <- (data.perc[, year_cols] / as.numeric(data.perc$DBH)) * 100
# Reality check to see that percentages add up to 100
all.equal(current = rowSums(data.perc[, year_cols], na.rm = TRUE),
          target = rep(100, length = nrow(data.perc)))






#########LIDAR PROCESSING################
##This is a 'meta-package' that lists all of the R packages for lidar and forest analysis. 
#devtools::install_github("atkinsjeff/ForestAnalysisInR") 
#library(ForestAnalysisInR)
#library(shiny)
#launchRFA()


##install TreeLS package that is not on CRAN
#remotes::install_github('tiagodc/TreeLS')

##install latest version of lidr tools (bug fixes) from github. Latest version is not on CRAN yet. 
#remotes::install_github('Jean-Romain/lidr')


library(raster)
library(lidR)
library(TreeLS)
library(rgdal) 
library(sf) 
library(rgeos) 
library(sp)
#lidR_parallelism

set_lidr_threads(30)
get_lidr_threads()

##Bring las point cloud into Rstudio
points = lidR::readLAS("D:/projects/Babst_TLidar_Tree_Rings/shared/pointcloud/Upper_Baldy_SantaRita/TLS_SantaRita_UBplot_pointcloud_041922_UTM_CC.las")

##Bring FARO Scene georeferenced point cloud (ASCII) into R. There are strange differences between the ascii and las point cloud. 
##These differences are impacted the analysis
#points = readTLS("UB_pointcloud_xyz_ascii.xyz")



##Voxelize to reduce point cloud density to 1 pnt/3 cubic cm
points_voxel = voxelize_points(points, 0.03)



##Delete the points on periphery of the cloud that are outside of the 'plot'
##This clipping procedure will be based on user defined plot center coordinate (GPS) and radius of plot
#The plot center coordinates and radii are currently found in the file 'Site_overview_SantaRitaMts_October2021.xls

center_coord_easting_x = 513637.25

center_coord_northin_y = 3507116.26

points_voxel_clip = clip_circle(points_voxel, center_coord_easting_x, center_coord_northin_y, 30)


##Separate trees from ground using lidr
points_classify = classify_ground(points_voxel_clip, algorithm = csf(sloop_smooth = TRUE, class_threshold = 0.1, cloth_resolution =  0.5, rigidness = 2))

##Make a point cloud file for just ground points
ground_points = filter_poi(points_classify, Classification == 2)

#Make a point cloud file for just non-ground (trees) points
tree_points = filter_poi(points_classify, Classification != 2)

##Create a grided digital terrain model from ground points
DTM = grid_terrain(ground_points, res = 0.5, algorithm = knnidw(k = 10, p = 2))

##Calculate vertical distance of each tree point above the DTM. This creates a vegetation height point cloud.
normalized_height = normalize_height(las = tree_points, algorithm = DTM)

##Remove points that are below 0
AGL_clean = filter_poi(normalized_height, Z > 0.5)

##Export normalized point cloud to .las
#writeLAS(AGL_clean, "F:\\Black_mypassport_1.8TB\\Smith_vanLeeuwen_work\\Babst_T-lidar\\tlidar_workspace\\UB_AGL.las")


##Create canopy height model raster. I don't think this is totally necessary. 
#CHM = rasterize_canopy(AGL_clean, res = 0.05)




##Identify stems with 'TreeLS'

##Identify tree occurrences for a normalized point cloud. It uses a Hough Transform which is a circle search
##Merge is a distance between stems. If two occurrences are less than the specified distance, they will be merged. 
##I don't completely understand how this algorithm works and what the parameters mean yet. It does seem to work, though. 

#tree_map2 = treeMap(AGL_clean, method = map.hough(min_h = 1, max_h = 3, h_step = 0.5, pixel_size = 0.025, 
# max_d = 0.85, min_density = 0.1, min_votes = 3), merge = 0.2, positions_only = FALSE)

#The eigen method seems to find more tree occurrences than the hough method
tree_map_eigen = treeMap(AGL_clean, method = map.eigen.knn(max_curvature = 0.1,
                                                           max_verticality = 10,
                                                           max_mean_dist = 0.1,
                                                           max_d = 0.5,
                                                           min_h = 1.5,
                                                           max_h = 3), merge = 0.2, positions_only = FALSE)


##Plot the voxelized point cloud with the identified tree stems
#AGL_clean_topped = filter_poi(AGL_clean, Z < 5)
x = plot(AGL_clean, color = "Z")
add_treeMap(x, tree_map_eigen, color = 'yellow', size = 3)



##Count the number of unique trees
unique_trees = unique(tree_map_eigen$TreeID)
length(unique_trees)

##Create 2D point layer showing the locations of the tree stems
#xymap = treeMap.positions(tree_map_eigen)

## Export 2D tree stem location points in a CSV
#write.csv(xymap,"F:\\Black_mypassport_1.8TB\\Smith_vanLeeuwen_work\\Babst_T-lidar\\tlidar_workspace\\UB_stems_eigen_29april2022.csv", row.names = TRUE)



##Identify all of the points that belong to each tree occurrence found in 'treeMap'
tree_points = treePoints(AGL_clean, tree_map_eigen, trp.voronoi())

##Plot AGL_clean with tree points and tree IDs
x = plot(AGL_clean, color = "Z")
add_treePoints(x, tree_points, size=4)
add_treeIDs(x, tree_points, color='yellow', cex=2)



##Identify stem points in each segmented tree
stem_points_eigen = stemPoints(tree_points, method = stm.eigen.knn(h_step = 0.5,
                                                                   max_curvature = 0.1,
                                                                   max_verticality = 10,
                                                                   voxel_spacing = 0.025,
                                                                   max_d = 0.5,
                                                                   votes_weight = 0.2,
                                                                   v3d = FALSE))

stem_points_hough = stemPoints(tree_points, method = stm.hough(h_step = 0.5,
                                                               max_d = 0.5,
                                                               h_base = c(1, 2.5),
                                                               pixel_size = 0.025,
                                                               min_density = 0.1,
                                                               min_votes = 3))




#plot AGL_clean point cloud with stem points
x = plot(AGL_clean, color = "Z")
add_stemPoints(x, stem_points_hough, color='red', size=8)

#plot only stem points
stem_points_hough_isolated = filter_poi(stem_points_hough, Stem == TRUE)
plot(stem_points_hough_isolated , color = "Stem")


#plot AGL_clean point cloud with stem points
x = plot(AGL_clean, color = "Z")
add_stemPoints(x, stem_points_eigen, color='red', size=8)

#plot only stem points
stem_points_eigen_isolated = filter_poi(stem_points_eigen, Stem == TRUE)
plot(stem_points_eigen_isolated , color = "Stem")

#writeLAS(tree_map2, "F:\\Black_mypassport_1.8TB\\Smith_vanLeeuwen_work\\Babst_T-lidar\\tlidar_workspace\\UB_tree_map.las")



##fit many cylinders to stem points
stem_seg = stemSegmentation(stem_points_hough_isolated, sgt.ransac.circle(tol = 0.1, n = 10, conf = 0.95, inliers = 0.9))

##plot stem segments
#x = plot(AGL_clean)
add_stemSegments(x, stem_seg, color='blue')
add_treeIDs(x, tree_points, color='yellow', cex=1)

tlsPlot(stem_seg)

##plot a specific tree segment
tlsPlot(stem_seg, fast = FALSE, segment = 5)



inv = tlsInventory(stem_points, dh = 1.3, dw = 0.5, hp = 1, d_method = shapeFit(shape = "circle", 
                                                                                algorithm = "ransac", n = 15, n_best = 20))


plot(stem_seg)


shape_fit = shapeFit(stem_segment = NULL,
                     shape = "circle",
                     algorithm = "ransac",
                     n = 10,
                     conf = 0.95,
                     inliers = 0.9,
                     n_best = 10,
                     z_dev = 30)



writeLAS(stem___points, "F:\\Black_mypassport_1.8TB\\Smith_vanLeeuwen_work\\Babst_T-lidar\\tlidar_workspace\\stem___points.las")






##Cindy, can you match the ground reference tree locations with the treeMap positions using a spatial join? 

######TREERING_GROUND_LIDAR Spatial Join###############
UB_plot_ground <- read.csv("D:/projects/Babst_TLidar_Tree_Rings/CLNorton/TLS_SantaRita_Plot_UB/UB_UTM.csv")
UB_plot_lidar <- read.csv("D:/projects/Babst_TLidar_Tree_Rings/UB_stem_locations.csv")


###tree_ring_ground_spatialpoints 1
xy_ground <- data.frame(treeID = UB_plot_ground$plotID, X = UB_plot_ground$lon, Y = UB_plot_ground$lat) ##CONVERTLATTOLONG
coordinates(xy_ground) <- ~ X + Y
proj4string(xy_ground) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  ## for example
tree.points.ground <- spTransform(xy_ground, CRS("+proj=utm +zone=12 +ellps=WGS84")) %>% st_as_sf()





###tree_ring_ground_spatialpoints 2
xy_ground <- data.frame(treeID = UB_plot_ground$plotID, X = UB_plot_ground$X, Y = UB_plot_ground$Y)
coordinates(xy_ground) <- c("X", "Y")
proj4string(xy_ground) <- CRS("+proj=utm +zone=12 +ellps=WGS84")  ## for example
tree.points.ground <- spTransform(xy_ground, CRS("+proj=utm +zone=12 +ellps=WGS84"))%>% st_as_sf()


##lidarspatiapoints
xy_lidar <- data.frame(lidar_treeID = UB_plot_lidar$TreeID, X = UB_plot_lidar$X, Y = UB_plot_lidar$Y)
coordinates(xy_lidar) <- c("X", "Y")
proj4string(xy_lidar) <- CRS("+proj=utm +zone=12 +ellps=WGS84")  ## for example
tree.points.lidar <- spTransform(xy_lidar, CRS("+proj=utm +zone=12 +ellps=WGS84")) %>% st_as_sf()




st_crs(tree.points.ground) = st_crs(tree.points.lidar)
tree.points.lidar_buffer = st_buffer(tree.points.lidar, 1.7)
tree_intersection = st_intersection(tree.points.lidar_buffer, tree.points.ground)




















##LIDAR
#Bring point cloud into Rstudio
points_georeferenced <- lidR::readLAS("D:/projects/Babst_TLidar_Tree_Rings/CLNorton/TLS_SantaRita_Plot_UB/TLS_SantaRita_UBplot_pointcloud_041922_aligned_CC.las",  filter = "-drop_class 7")
points <- lidR::readLAS("D:/projects/Babst_TLidar_Tree_Rings/CLNorton/TLS_SantaRita_Plot_UB/TLS_SantaRita_UBplot_pointcloud_021622_17error.las")

#metadata <-readxl::read_xlsx("D://projects//Babst_TLidar_Tree_Rings//CLNorton//TLS_SantaRIta_Plot_UB//Metadata_SantaRitaMts_Winter2021_2022.xlsx", sheet = 1)
#tree_metrics = function(point_cloud,tree_ring, ground_data) {
#identify point cloud center
xcenter = 513637.25
ycenter = 3507116.26
#clip raster based on center and radius
clip_points_georeferenced <- clip_circle(points_georeferenced, radius = 30, xcenter = xcenter, ycenter =	ycenter)

#writeLAS(clip_points_georeferenced, file = "clip_points_georeferenced_042022.las")
#plot(clip_points_georeferenced)
#Identify ground points using a cloth simulation filter algorithm
points_classified <- classify_ground(clip_points_georeferenced, algorithm = csf(sloop_smooth = FALSE, class_threshold = 0.4, cloth_resolution =  0.25, rigidness = 2))

#Separate the classified point cloud into a ground point cloud and canopy point cloud
ground_points <- filter_poi(points_classified, Classification == 2)
canopy_points <- filter_poi(points_classified, Classification == 0)
#writeLAS(canopy_points, file = "canopy_points_042022.las")
#Create a grided digital terrain model
DTM = grid_terrain(ground_points, res = 0.25, algorithm = knnidw(k = 10, p = 2), na.omit())

#writeRaster(DTM, file = "DTM_042022.tif")
#Calculate vegetation height above the ground
AGL <- normalize_height(canopy_points, DTM, na.rm = TRUE, copy=TRUE)

#Remove points that are below 0
AGL_clean <- filter_poi(AGL,Z > 0 | Z <35)

#Create canopy height model(raster) from the normalized canopy points
CHM <- grid_canopy(AGL_clean, res = 0.25, algorithm = p2r(subcircle = 0, na.fill = NULL))
#Voxelize the canopy point cloud. This thins the point cloud while retaining structure. It speeds processing commands later in the code. 
#Each point represents a 25 cm cube
canopy_voxel <- voxelize_points(AGL_clean, res = 0.25)
#Identify individual trees. Tree tops are identified by a user defined moving window looking for local maximums across CHM

#lidr segmentation and find trees
treetops_lidR <- find_trees(canopy_voxel, lmf(ws=8, shape = "circular"))
#treemap_TreeLiDAR = treeMap(AGL_clean, method = map.hough(min_h = 1, max_h = 3, h_step = 0.5, pixel_size = 0.025, 
#                                                 max_d = 0.85, min_density = 0.1, min_votes = 3), merge = 0.2, positions_only = FALSE)
#unique_trees = unique(tree_map$TreeID)
#xymap = treeMap.positions(tree_map) 

##lidar_spatialpoints
#xy_lidar <- data.frame(lidar_treeID = xymap$TreeID, X = xymap$X, Y = xymap$Y)
#coordinates(xy_lidar) <- c("X", "Y")
#proj4string(xy_lidar) <- CRS("+proj=utm +zone=12 +ellps=WGS84")  ## for example
#treemap_TreeLiDAR_SPDF <- spTransform(xy_lidar, CRS("+proj=utm +zone=12 +ellps=WGS84"))

TLS_lidR <- segment_trees(AGL_clean, dalponte2016(chm = CHM, treetops = treetops_lidR, th_tree = 8, th_seed = 0.20,
                                                  th_cr = 0.20, max_cr = 35, ID = "treeID"))

#plot(TLS_lidR)
##Plot the voxelized point cloud with the identified tree stems
#x = plot(AGL_clean, color = "Z")
#add_treeMap(x, tree_map2, color = 'yellow', size = 2)
##Count the number of unique trees
#length(unique_trees)
##Create 2D point layer showing the locations of the tree stems
xymap = treeMap.positions(tree_map) 
##Create canopy height model raster. I don't think this is totally necessary. 
CHM = rasterize_canopy(AGL_clean, res = 0.05)


#Calculate a convex hull and draw a polygon for each tree
metric <- tree_metrics(TLS_lidR, .stdtreemetrics)
crown_Hull  = delineate_crowns(TLS_lidR)
crown_Hull@data = dplyr::left_join(crown_Hull@data, metric@data)

plot(crown_Hull)

####FUNCTION 5#####
##JGillan
####LIDAR PROCESSING######
set_lidr_threads(30)
get_lidr_threads()

####PREPARE POINT CLOUD####
##Delete the points on periphery of the cloud that are outside of the 'plot'
##This clipping procedure will be based on user defined plot center coordinate (GPS) and radius of plot
#The plot center coordinates and radii are currently found in the file 'Site_overview_SantaRitaMts_October2021.xls
center_coord_easting_x = 513637.25
center_coord_northin_y = 3507116.26

points_voxel = voxelize_points(lidar, 0.03)
points_voxel_clip = clip_circle(points_voxel, center_coord_easting_x, center_coord_northin_y, 30)
##Separate trees from ground using lidr
points_classify = classify_ground(points_voxel_clip, algorithm = csf(sloop_smooth = TRUE, class_threshold = 0.1, cloth_resolution =  0.5, rigidness = 2))
##Make a point cloud file for just ground points
ground_points = filter_poi(points_classify, Classification == 2)
#Make a point cloud file for just non-ground (trees) points
tree_points = filter_poi(points_classify, Classification != 2)
##Create a grided digital terrain model from ground points
DTM = grid_terrain(ground_points, res = 0.5, algorithm = knnidw(k = 10, p = 2))
##Calculate vertical distance of each tree point above the DTM. This creates a vegetation height point cloud.
normalized_height = normalize_height(las = tree_points, algorithm = DTM)
##Remove points that are below 0
AGL_clean = filter_poi(normalized_height, Z > 0.5)


##Export normalized point cloud to .las
#writeLAS(AGL_clean, "F:\\Black_mypassport_1.8TB\\Smith_vanLeeuwen_work\\Babst_T-lidar\\tlidar_workspace\\UB_AGL.las")


###ITD###
##Identify stems with 'TreeLS'
##Identify tree occurrences for a normalized point cloud. It uses a Hough Transform which is a circle search
##Merge is a distance between stems. If two occurrences are less than the specified distance, they will be merged. 
##I don't completely understand how this algorithm works and what the parameters mean yet. It does seem to work, though. 

tree_map = treeMap(AGL_clean, method = map.hough(min_h = 1, max_h = 3, h_step = 0.5, pixel_size = 0.025, 
                                                 max_d = 0.85, min_density = 0.1, min_votes = 3), merge = 0.2, positions_only = FALSE)


#plot(TLS_lidR)
##Plot the voxelized point cloud with the identified tree stems
#x = plot(AGL_clean, color = "Z")
#add_treeMap(x, tree_map2, color = 'yellow', size = 2)
##Count the number of unique trees
unique_trees = unique(tree_map$TreeID)
#length(unique_trees)
##Create 2D point layer showing the locations of the tree stems
xymap = treeMap.positions(tree_map) 
##Create canopy height model raster. I don't think this is totally necessary. 
CHM = rasterize_canopy(AGL_clean, res = 0.05)
#Tip: to grow the regions bigger, make th_seed and th_cr have small values (e.g. < 0.25)
TLS_lidR <- segment_trees(AGL_clean, dalponte2016(chm = CHM, treetops = tree_map, th_tree = 8, th_seed = 0.20,
                                                  th_cr = 0.20, max_cr = 35, ID = "treeID"))


######LIDAR TREE METRICS and lidar can be merged with TLS and ground from here
#Calculate a convex hull and draw a polygon for each tree
metric <- tree_metrics(TLS_lidR, .stdtreemetrics)
crown_Hull  = delineate_crowns(TLS_lidR)
crown_Hull@data = dplyr::left_join(crown_Hull@data, metric@data)


#calculate crown diameter
crown_Hull[["crown_Diameter"]] = sqrt(crown_Hull[["convhull_area"]]/ pi) * 2 #crown diameter

#calculate tree volume
#tree volume base number of voxels multiplied by the volume of voxels assigned
crown_Hull[["tree_volume_cm3"]] = (crown_Hull[["npoints"]]*25)  #vokume cm3
crown_Hull[["tree_Volume_m3"]] = (crown_Hull[["npoints"]]*25)/1000000  #volume m3


#crown_Hull summary based on area and height
sp_summarise(crown_Hull, variables = c("convhull_area", "Z")) 

#put in data frame and rename columns
tree_dataframe = as.data.frame(crown_Hull) %>% rename(Height = Z, crown_Area = convhull_area)
tree_dataframe

##EXPORT
#Export the convex hull polygons into a shapefile 
writeOGR(crown_Hull, dsn = "D://projects//Babst_TLidar_Tree_Rings//CLNorton//TLS_SantaRIta_Plot_UB//Output//030222", layer = "crown_Hull_UB", driver = "ESRI Shapefile", overwrite_layer = TRUE)

#Export the treetop points into a shapefile. This can be used to understand how trees were identified. Not totally necessary in the workflow. 
writeOGR(treetops, dsn = "D://projects//Babst_TLidar_Tree_Rings//CLNorton//TLS_SantaRIta_Plot_UB//Output//030222", layer = "treetops_UB", driver = "ESRI Shapefile", overwrite_layer = TRUE)

#Export the Canopy Height Model to a .tif
writeRaster(CHM, "D://projects//Babst_TLidar_Tree_Rings//CLNorton//TLS_SantaRIta_Plot_UB//Output//030222CHM_UB", overwrite = TRUE)

write.csv(tree_dataframe, file = "UB_tree_TLS_dataframe.csv")






#ground validation parsing draft
ground_tree_volume <- function(metadata) {
  volume <- metadata %>% select(c(DBH,Height)) %>%
    transform(DBH_in = DBH/2.54)%>%
    transform(Height_feet = Height*3.28084)%>%
    transform(board_16feet = Height_feet/16)
  transform(volume = pi*DBH_in*Height)
  
  return(volume)
  
}




#DBH SCRIPT DRAFT - segmentation - 16 foot foot slices
##try filtering points after tree dileanation segmentation
bh_points <- filter_poi(TLS_dalponte, Z < 1.5)

#Each point represents a 25 cm cube
bh_voxel <- voxelize_points(bh_points, res = 0.25)

bh_tops <- find_trees(bh_voxel, lmf(ws=0.1, shape = "circular"))

bh_CHM <- grid_canopy(bh_points, res = 0.25, algorithm = p2r(subcircle = 0, na.fill = NULL))

#Tip: to grow the regions bigger, make th_seed and th_cr have small values (e.g. < 0.25)
TLS_dalponte <- segment_trees(bh_voxel, dalponte2016(chm = bh_CHM, treetops = bh_tops, th_tree = 1, th_seed = 0.1,
                                                     th_cr = 0.1, max_cr = 35, ID = "treeID"))
#Calculate a convex hull and draw a polygon for each tree
metric <- tree_metrics(bh_points, .stdtreemetrics)
metric$convhull_area

bhHull  = delineate_crowns(bh_points)
bhHull@data = dplyr::left_join(bhHull@data, metric@data)
plot(bhHull)
plot(crown_Hull)
#calculate bh diameter
bhHull[["bh_Diameter"]] = sqrt(bhHull[["convhull_area"]]/ pi) * 2 #bh diameter

#bhHull summary based on area and height
sp_summarise(bhHull, variables = c("convhull_area", "Z")) 

#put in data frame and rename columns
dbh_dataframe = as.data.frame(bhHull) %>% rename(Height = Z, bhArea = convhull_area, )

dbh_dataframe
tree_dataframe










##TREERING
Hainich_Treerings <- read.rwl("D:/projects/Babst_TLidar_Tree_Rings/CLNorton/Europe TLS Data/TLS Data/Hainich_Treerings/Good_FASY.rwl") %>% as.data.frame()
mean_width_yr <- Hainich_Treerings %>% rowMeans(na.rm = TRUE) %>% as.data.frame() %>% rownames_to_column()
plot(mean_width_yr,main="Hainich",
     xlab="Year",
     ylab="growth increments",
     cex.axis=1.5,
     font.main=2, font.lab=2, font.sub=2,
     cex.main=3, cex.lab=1.7, cex.sub=1.2)
lines(mean_width_yr$rowname,mean_width_yr$., col="red", lwd=2)
abline(h=mean(mean_width_yr$.), col="black", lwd=3, lty=2)

p <- ggplot(data = mean_width_yr,                # specify dataset
            aes(x = rowname, y = .)) + # Age on x-, Raven on y-axis
  geom_point(pch = 1)                  # plot points (pch = 1: circles, type '?pch' for other options)
p

p <- p +
  xlab("Age (years)") +            # add label for x axis
  ylab("Raven's Matrices (score)") # add label for y axix
p
q1 <- p +
  geom_smooth() # add scatterplot smoother
q1





writeOGR(treetops_lidR, ".", "TLS_lidR_ITD", 
         driver = "ESRI Shapefile") 
writeOGR(crown_Hull, dsn = "D://projects//Babst_TLidar_Tree_Rings//CLNorton//TLS_SantaRIta_Plot_UB//Output//030222", layer = "crown_Hull_UB_lidR", driver = "ESRI Shapefile", overwrite_layer = TRUE)



#watershed ForestTools segmentation
lin = function(x){x * 0.05 + 0.6}
#filter low-lying underbrush and spurious treetops
treetops_ForestTools = vwf(CHM = CHM, winFun = lin, minHeight = 8)
TLS_ForestTools = mcws(treetops = treetops_lidR, CHM = CHM,format = "polygons", minHeight = 8, verbose = FALSE)

#plot(TLS_ForestTools, col = sample(rainbow(50), length(unique(TLS_ForestTools[])), replace = TRUE), legend = FALSE, xlab = "", ylab = "", xaxt='n', yaxt = 'n')
plot(TLS_ForestTools, border = "blue", lwd = 0.5, add = TRUE)

writeRaster(CHM, file = "CHM_UB.tif")


writeOGR(treetops_ForestTools, ".", "TLS_ForestTools_ITD", 
         driver = "ESRI Shapefile") 
writeOGR(TLS_ForestTools, ".", "TLS_ForestTools_crown", 
         driver = "ESRI Shapefile") 




##rLiDAR
# Smoothing CHM
schm<-CHMsmoothing(CHM, "mean", 5)
# Setting the fws:
fws<-8 # dimention 5x5
# Setting the specified height above ground for detectionbreak
minht<-10
# Getting the individual tree detection list
treeList<-FindTreesCHM(schm, fws, minht)
plot(treeList)
summary(treeList)
# Plotting the individual tree location on the CHM
plot(CHM, main="LiDAR-derived CHM")
XY<-SpatialPoints(treeList[,1:3]) # Spatial points
treetops_rLiDAR <- SpatialPointsDataFrame(XY, data.frame(row.names=row.names(XY),
                                                         ID=1:length(XY)))

TLS_rLiDAR <- ForestCAS(chm = schm, loc = as.data.frame(XY))
plot(TLS_rLiDAR)

boundaryTrees<-TLS_rLiDAR[[1]]
# Plotting the individual tree canopy boundary over the CHM
plot(chm, main="LiDAR-derived CHM")
# adding tree canopy boundary
plot(boundaryTrees, add = TRUE, border = 'red', bg = 'transparent')




writeOGR(treetops_rLiDAR, ".", "TLS_rLiDAR_ITD", 
         driver = "ESRI Shapefile") 
writeOGR(boundaryTrees, ".", "TLS_rLiDAR_crown", 
         driver = "ESRI Shapefile") 



#Export the convex hull polygons into a shapefile 
writeOGR(hulls4, dsn = "D:/projects/Babst_TLidar_Tree_Rings/lidar_tests/Output", layer = "TLS_trees_UBplot", driver = "ESRI Shapefile", overwrite_layer = TRUE)

#Export the treetop points into a shapefile. This can be used to understand how trees were identified. Not totally necessary in the workflow. 
writeOGR(treetops, dsn = "D:/projects/Babst_TLidar_Tree_Rings/lidar_tests/Output", layer = "TLS_treetops_plot", driver = "ESRI Shapefile", overwrite_layer = TRUE)

#Export the Canopy Height Model to a .tif
writeRaster(CHM, "D:/projects/Babst_TLidar_Tree_Rings/lidar_tests/Output/TLS_CHM_UBplot.tif", overwrite = TRUE)







#FORESTTOOLS
#install.packages("ForestTools")
#crowndiameter and ITD

# Attach the 'ForestTools' and 'raster' libraries
library(ForestTools)
library(raster)

# Load sample canopy height model
data("CHM")


# Remove plot margins (optional)
par(mar = rep(0.5, 4))

# Plot CHM (extra optional arguments remove labels and tick marks from the plot)
plot(CHM, xlab = "", ylab = "", xaxt='n', yaxt = 'n')

lin = function(x){x * 0.05 + 0.6}

#filter low-lying underbrush and spurious treetops
ttops = vwf(CHM = CHM, winFun = lin, minHeight = 0.5)
windows()
plot.new()
# Add dominant treetops to the plot
plot(ttops, col = "blue", pch = 20, cex = 0.5, add = TRUE)


mean(ttops$height)


# Create crown map
crowns = mcws(treetops = ttops, CHM = CHM, minHeight = 1, verbose = FALSE)

# Plot crowns
plot(crowns, col = sample(rainbow(50), length(unique(crowns[])), replace = TRUE), legend = FALSE, xlab = "", ylab = "", xaxt='n', yaxt = 'n')

crowns$layer

# Create polygon crown map
crownsPoly = mcws(treetops = ttops, CHM = CHM, format = "polygons", minHeight = 1, verbose = FALSE)

crownsPoly$height
# Add crown outlines to the plot
plot.new()
plot(crownsPoly, border = "blue", lwd = 0.01, add = TRUE)
crownsPoly$

# Compute average crown diameter
crownsPoly[["crownDiameter"]] = sqrt(crownsPoly[["crownArea"]]/ pi) * 2

# Mean crown diameter
mean(crownsPoly$crownDiameter)
sp_summarise(ttops)
sp_summarise(crownsPoly, variables = c("crownArea", "height"))
data("Blocks")

# Compute tree count and height statistics for cut blocks
blockStats = sp_summarise(ttops, areas = Blocks, variables = "height")

# Plot CHM
plot(CHM, xlab = "", ylab = "", xaxt='n', yaxt = 'n')

# Add block outlines to the plot
plot(Blocks, add = TRUE, border =  "darkmagenta", lwd = 2)

# Add tree counts to the plot
library(rgeos)
text(gCentroid(Blocks, byid = TRUE), blockStats[["TreeCount"]], col = "darkmagenta", font = 2)

blockStats@data





##HENDRIX



# Clear workspace:
rm(list=ls())

library(lidR)
library(feather)
library(data.table)
#terminal output coloring
library(crayon)
error = red $ bold
warn = magenta $ underline
note = cyan
#cat(error("Error: subscript out of bounds!/n"))
#cat(warn("Warning: shorter argument was recycled./n"))
#cat(note("Note: no such directory./n"))



#then change working directory for writing output:
setwd("D:/projects/Babst_TLidar_Tree_Rings/lidar_tests/Output")
las_classified <- classify_ground(points, algorithm = csf(sloop_smooth = TRUE, class_threshold = 0.55, cloth_resolution =  0.2, rigidness = 2))

#make dtm:
dtm <- grid_terrain(las_classified, res = .1, algorithm = knnidw(k = 8, p = 2), na.rm = TRUE)
dtmRaster <- raster::as.raster(dtm)
raster::plot(dtmRaster)

#lasRaster = raster::as.raster(las)

lasnorm = normalize_height(las_classified, dtm)

# compute a canopy image
chm = grid_canopy(lasnorm, res = 0.1, algorithm = p2r(subcircle = 0, na.fill = NULL)) 
chm = raster::as.raster(chm)
dev.new()
raster::plot(chm)

dev.off()

# smoothing post-process (e.g. 2x mean)
kernel = matrix(1,5,5)
schm = focal(chm, w = kernel, fun = mean)
#schm = raster::focal(chm, w = kernel, fun = mean)

dev.new()
raster::plot(chm, col = height.colors(50)) # check the image
quartz.save("LS_SCHM.png")

# save smoothed canopy height model as tif
raster::writeRaster(schm, "TLS_Outliers_Removed_Smoothed_CHM_no_edge_stretch.tif", format = "GTiff", overwrite = TRUE)

# tree segmentation
# 'th' Numeric. Number value below which a pixel cannot be a crown.
#Default 2
crowns  segment_trees(lasnorm, watershed(chm = schm, treetops = treetops, th_tree = 1.5, th_seed = 0.20,
                                                  th_cr = 0.20, max_cr = 35, ID = "treeID"))

# Plotting point cloud of trees only:
# without rendering points that are not assigned to a tree
tree = lasfilter(lasnorm, !is.na(treeID))
# this would be a good plot to play on rotate using ImageMagick:
plot(tree, color = "treeID", colorPalette = pastel.colors(100))

#save tree point cloud (clustered point cloud):
writeLAS(tree, "A-lidar_Clustered_By_Watershed_Segmentation.las")
#write.csv(tree@data, "all20TilesGroundClassified_and_Clustered_By_Watershed_Segmentation.csv")
write_feather(tree@data, "A-lidar_Clustered_By_Watershed_Segmentation.feather")

# Plotting raster with delineated crowns:
library(raster)
contour = rasterToPolygons(crowns, dissolve = TRUE)

dev.new()
plot(schm, col = height.colors(50))
plot(contour, add = T)
quartz.save("A-lidar Segmented SCHM - MCC-Lidar Classing & KNN-IDW Rasterization.png")
#dev.off()
