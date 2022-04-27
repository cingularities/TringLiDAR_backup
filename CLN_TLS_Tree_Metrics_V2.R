library(dplR)
library(tidyr)
library(tidyverse)

setwd("D:/projects/Babst_TLidar_Tree_Rings/CLN_TLS")



############################################ Tree Ring and Ground Fusion ####################################################

#select the path for folder where you have the file to load
ground <- "D:/projects/Babst_TLidar_Tree_Rings/CLN_TLS/Europe TLS Data/Serial_sampling/metadata_hainich_Mar2019.xlsx"
tree<- "D:/projects/Babst_TLidar_Tree_Rings/CLN_TLS/Europe TLS Data/TLS Data/Hainich_Treerings/Good_FASY.rwl"
#lidar = lidR::readLAS("D:/projects/Babst_TLidar_Tree_Rings/TLS_SantaRita_UBplot_pointcloud_041922_aligned_CC.las")



###FUNCTION 1 - tree ring data parse and formatting###
#raw tree ring data parsing and mean tree cores
treering_parse <- function(tree) {
  rwl<-dplR::read.rwl(tree) #read file path of rwl
  detrend.rwi <- dplR::detrend(rwl = rwl, method = "Spline") %>% #detrend rwl
    t() %>% #transpose years as columns and row as each core sample 
    as.data.frame() %>% #as data frame
    dplyr::rename_with( ~ paste0("year_", .x)) %>% #adding "year_" to each year column name for unique identifiers
    rownames_to_column()%>% #converts sample row name into a column to be parsed
    tidyr::separate(rowname,
                    into = c("site", "coreID"),
                    sep = "(?<=[A-ZA-Z])(?=[0-9])") %>% ## separates column of site and treeID
    tidyr::separate(coreID,
                    into = c("treeID", "core"),
                    sep = "(?<=[0-9])(?=[A-ZA-Za-z])")  ##additional splot of CoreID to separate core sample from tree
  return(detrend.rwi)}
treering_function_results <- treering_parse(tree)

#write.csv(data.transpose, file = "datatranspose.csv")










##RAW GROUND DATA
coor <- c(51.079167, 10.453)
tree23 <- destPoint(coor,110, 76)
plotcenter <- destPoint(tree23,98, 28.1)
hainich_center<-c(51.08007, 10.45273)

###FUNCTION 2 - ground data azimuth and distance to coordinates###
ground_parse <- function(coor, ground) {
  
  dist_azim <- readxl::read_excel(ground) #read xlsx ground metadata file

  coor_list <- list() #creates empty list for new coordinates
  
  for (i in 1:nrow(dist_azim)) { #for loop to create a lat long for each azimuth and distance from plot center lat and lon
    stemAzimuth <- dist_azim$Azimuth[i] #grab azimuth
    stemDistance <- dist_azim$Distance[i] #grab distance
    
    new_coordinates <- geosphere::destPoint(coor, stemAzimuth,  stemDistance) %>% #estimates new lat lon coordinate
      as.data.frame() %>% #make as data frame
      add_column(treeID = dist_azim$TreeID[i]) %>% #add treeID column of that tree
      add_column(stemAzimuth = stemAzimuth) %>% #add Azimuth column of that tree
      add_column(stemDistance = stemDistance) %>% #add Distance column of that tree
      add_column(species = dist_azim$Species[i]) %>% #add Species column of that tree
      add_column(Height = dist_azim$Height[i])%>% #add Height column of that tree
      add_column(DBH = dist_azim$DBH[i]) #add DBH column of that tree
    
    
    coor_list[[i]] <- new_coordinates #put result of loop in list
    
  }
  
  coor_df<- coor_list %>% bind_rows() %>% na.omit() #remove NA rows
  
  
  return(coor_df)
}

ground_function_results <- ground_parse(hainich_center, ground) 














###FUNCTION 3 -  tree ring and ground data merge###
###TREE rings and ground merge#####
treering_ground <- function (treering_parse, ground_parse, plot) {
  
  treering_ground_join <- treering_parse %>% 
    subset(site == plot) %>% #subsets to input plot
    lapply(as.numeric)%>% #makes numeric
    as.data.frame() %>% #makes into df
    group_by(treeID) %>% #groups data by the treeID
    summarise(across(-c(site, core), mean, na.rm = TRUE)) %>% #a mean summary of multiple cores per tree
    rev()%>% #rev the entire data frame
    left_join(ground_parse, by = "treeID") #joins tree ring data with ground metadat using treeID
  
  return(treering_ground_join)
}

treering_ground_results <- treering_ground(treering_function_results, ground_function_results, "HF")
























###FUNCTION 4 - DBH and Percentage reconstruction###
###DBH RECONSTRUCITONS
DBH_recon <- function(treering_ground_join){
  
  treering_ground_join[is.na(treering_ground_join)] <- 0  #zero all NA
  year_cols <- which(substr(x = colnames(treering_ground_join), start = 1, stop = 4) == "year")  # Identify the columns that have data of interest which are year columns
  treering_ground_join$lp <- rowSums(x = treering_ground_join[, year_cols], na.rm = TRUE)  #row sums all the years into column lp, removing NA 
  lh <- treering_ground_join  # Make a copy of that data frame that we will use for differences
  lh[, year_cols] <- sapply(1:ncol(lh[, year_cols]), function(col){rowSums(lh[, year_cols][1:col])}) #convert years to cumulative sum of year columns
  diff <- lh # Make a copy of that data frame that we will use for difference
  diff[, year_cols] <- (diff$lp) - (diff[, year_cols]) #calculate difference
  fraction <- lh # Make a copy of that data frame that we will use for difference
  fraction[, year_cols] <- diff[, year_cols]/diff$lp #calculate difference
  diameter <- lh # Make a copy of that data frame that we will use for difference
  diameter[, year_cols] <- as.numeric(fraction$DBH)*fraction[, year_cols] #calculate difference
  
  # Reality check to see that percentages add up to 100
  
  return(diameter)
}



DBH_recon_results <- DBH_recon(treering_ground_results)



##PERCENTAGE DBH####
DBH_perc <- function(diameter) {
`year_cols <- which(substr(x = colnames(diameter), start = 1, stop = 4) == "year")# Identify the columns that have data of interest
 DBH_perc <- diameter # Make a copy of that data frame that we will use for percentages
 DBH_perc[, year_cols]<- (DBH_perc[, year_cols]/as.numeric(DBH_perc$DBH))* 100 # Do percentage calculations for only those columns with data of interest
return(DBH_perc)
}


DBH_perc_results <- DBH_perc(DBH_recon_results)














####FUNCTION 5#####
####LIDAR PROCESSING######
set_lidr_threads(30)
get_lidr_threads()


X <- X #X coordinate
Y <- Y #Y coordinate
lidar_voxel_clip = clip_circle(lidar_voxel, center_coordx, center_coordy, 25) #clip point cloud
lidar_classify = classify_ground(lidar_voxel_clip, 
                                 algorithm = csf(sloop_smooth = TRUE, 
                                                 class_threshold = 0.1, 
                                                 cloth_resolution =  0.5, 
                                                 rigidness = 2)) ##Separate trees from ground using lidr
ground_lidar = filter_poi(lidar_classify, Classification == 2) #Make a point cloud file for just ground lidar
tree_lidar = filter_poi(lidar_classify, Classification != 2) #Make a point cloud file for just non-ground (trees) lidar
DTM = grid_terrain(ground_lidar, res = 0.5, algorithm = knnidw(k = 10, p = 2)) #Create a grided digital terrain model from ground lidar
normalized_height = normalize_height(las = tree_lidar, algorithm = DTM) #Calculate vertical distance of each tree point above the DTM. This creates a vegetation height point cloud.
AGL_clean = filter_poi(normalized_height, Z > 0.5) #Remove lidar that are vertically below 0.5 m

#Identify stems with 'TreeLS'
#Identify tree occurrences for a normalized point cloud. It uses a Hough Transform which is a circle search
#Merge is a distance between stems. If two occurrences are less than the specified distance, they will be merged. 

tree_map = TreeLS::treeMap(AGL_clean, method = map.hough(min_h = 1, max_h = 3, h_step = 0.5, pixel_size = 0.025, 
                                                         max_d = 0.85, min_density = 0.1, min_votes = 3), merge = 0.2, positions_only = FALSE)


#Plot the voxelized point cloud with the identified tree stems
x = plot(AGL_clean, color = "Z")
add_treeMap(x, tree_map, color = 'yellow', size = 2)

#Count the number of unique trees
unique_trees = unique(tree_map$TreeID)
length(unique_trees)

#Create 2D point layer showing the locations of the trees
xy_lidar = treeMap.positions(tree_map)

##Cindy, can you match the ground reference tree locations with the treeMap positions using a spatial join? 

######TREERING_GROUND_LIDAR Spatial Join###############
UB_plot_ground <- read.csv("D:/projects/Babst_TLidar_Tree_Rings/CLN_TLS/TLS_SantaRita_Plot_UB/UB_UTM.csv")
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









