
setwd("D:/LANDSCAPE_GENOMICS/EggplantCWR/Eggplant")

library(factoextra)
library(FactoMineR)
library(ggplot2)

library(dplyr)
library(tidyr)
library(conflicted)
#samples
epafro_wr = read.csv('D:/LANDSCAPE_GENOMICS/EggplantCWR/EggplantAfro/afroep_all_env_variables2.csv')
epafro_wr
epafro_cwr = epafro_wr %>% dplyr::select(genotype, ind, sp_order, species ,lon, lat,country) %>% drop_na()

epafro_cwr
nrow(epafro_cwr)


#####ENVIRONMENTAL VARIABBLES####

library(raster)
library(sp)

# download the bioclim data:

bioclims=list.files( path='D:/DATA/Climate data/wc2.1_2.5m_bio', patter='.tif', full.names = TRUE)
clims = raster::stack(bioclims)

#extra wclim 2.5 resolution
#minimum temperature (°C)
tmin_files = list.files("D:/DATA/Climate data/wc2.1_2.5m_tmin/", pattern = ".tif", full.names = TRUE)
tmin_stack = raster::stack(tmin_files)
tmin = mean (tmin_stack)


#maximum temperature (°C)

tmax_files = list.files("D:/DATA/Climate data/wc2.1_2.5m_tmax/", pattern = ".tif", full.names = TRUE)
tmax_stack = raster::stack(tmax_files)
tmax = mean (tmin_stack)

#average temperature (°C)

tavg_files = list.files("D:/DATA/Climate data/wc2.1_2.5m_tavg/", pattern = ".tif", full.names = TRUE)
tavg_stack = raster::stack(tavg_files)
tavg = mean (tavg_stack)

#precipitation (mm) 

prec_files = list.files("D:/DATA/Climate data/wc2.1_2.5m_prec/", pattern = ".tif", full.names = TRUE)
prec_stack = raster::stack(prec_files)
prec = mean (prec_stack)


#solar radiation (kJ m-2 day-1)

srad_files = list.files("D:/DATA/Climate data/wc2.1_2.5m_srad/", pattern = ".tif", full.names = TRUE)
srad_stack = raster::stack(srad_files)
srad = mean (srad_stack)


#wind speed (m s-1)

wind_files = list.files("D:/DATA/Climate data/wc2.1_2.5m_wind/", pattern = ".tif", full.names = TRUE)
wind_stack = raster::stack(wind_files)
wind = mean (wind_stack)

#water vapor pressure (kPa)
vapr_files = list.files("D:/DATA/Climate data/wc2.1_2.5m_vapr/", pattern = ".tif", full.names = TRUE)
vapr_stack = raster::stack(vapr_files)
vapr = mean (vapr_stack)

#Elevation

Elev = raster("D:/DATA/Climate data/wc2.1_2.5m_elev/wc2.1_2.5m_elev.tif")

#Extra variables files
wcl_extra = raster::stack(tmin, tmax, srad, vapr, wind, Elev)
names(wcl_extra) = c ("tmin", "tmax", "srad", "vapr", "wind", 'elev')

wcl_extra

#wclim variables

wclm_all = raster::stack(clims, wcl_extra)

#soil data (ISRICsoilgrids) (https://files.isric.org/soilgrids/latest/data/)

#gdalwarp(t_srs="EPSG:4326", multi=TRUE, wm=200, 
#co=c("BIGTIFF=YES", "COMPRESS=DEFLATE", "TILED=TRUE"),
#tr=c(0.25,0.25), # Desired output resolution
#"/vsicurl?max_retry=3&retry_delay=1&list_dir=no&url=https://files.isric.org/soilgrids/latest/data/phh2o/phh2o_15-30cm_mean.vrt", # Input VRT
#"phh2o.tif") # Output file)
#phh2o_15_30cm
ph = raster("D:/LANDSCAPE_GENOMICS/Gap analysis/soil data/phh2o_15-30cm_SoilGrids.tif")

#nitrogen_15_30cm
nitro = raster("D:/LANDSCAPE_GENOMICS/Gap analysis/soil data/nitrogen_15-30cm_SoilGrids.tif")
plot(nitro)

#soil organic carbon_15_30cm
soc = raster("D:/LANDSCAPE_GENOMICS/Gap analysis/soil data/soc_15-30cm_SoilGrids.tif")

#cation exchange capacity_15-30cm
cec = raster("D:/LANDSCAPE_GENOMICS/Gap analysis/soil data/cec_15-30cm_SoilGrids.tif")

#clay
clay = raster("D:/LANDSCAPE_GENOMICS/Gap analysis/soil data/clay_15-30cm_SoilGrids.tif")
#sand
sand = raster("D:/LANDSCAPE_GENOMICS/Gap analysis/soil data/sand_15-30cm_SoilGrids.tif")
#silt
silt = raster("D:/LANDSCAPE_GENOMICS/Gap analysis/soil data/silt_15-30cm_SoilGrids.tif")
#ocd
ocd = raster("D:/LANDSCAPE_GENOMICS/Gap analysis/soil data/ocd_15-30cm_SoilGrids.tif")
#ocs
ocs = raster("D:/LANDSCAPE_GENOMICS/Gap analysis/soil data/ocs_0-30cm_SoilGrids.tif")

#stack soil files

sfiles = c(ph, nitro, soc, cec, clay, sand, silt, ocd, ocs)
srasters =  raster::stack(sfiles)
srasters
#plot(srasters)

#stack soil and bio rasters to make a predictors raster

srasters2 = raster::resample(srasters, clims) # change soil rasters to same resolution as bio
res(srasters2)
res(wclm_all)

wclm_all2 = raster::crop(wclm_all, extent(srasters2)) # make wclm and soil raster to same extent


env = raster::stack(wclm_all2, srasters2)
names(env)

#coordinates as spatialpoints

epafro_coord = epafro_cwr %>% dplyr::select(lon, lat)
#coordinates(epafro_coord) = c('lon', 'lat')
epafro_coord = SpatialPoints(coords = epafro_coord, proj4string=CRS("+proj=longlat +ellips=WGS84"))
epafro_coord 

library(rnaturalearth)
library(sp)

#world countries
sp::plot(ne_countries())
points(epafro_coord, col= 'red')

#get the environmental variables at the sample coordinates
epafroenv_var = raster::extract(env, epafro_coord)


epafro_env = cbind(epafro_cwr, epafroenv_var)
head(epafro_env)
length(which(is.na(epafro_env))) #check for NAs

epafro_env1 = epafro_env
epafro_env1
#colnames(ep_climates1)= c("genotype","VINO","taxon","lon","lat","country","c_code","continent",
#"AMT_bio1","MTWaQ_bio10","MTCoQ_bio11", "AP_bio12","PWeM_bio13",
#"PDrM_bio14" ,"PSe_bio15","PWeQ_bio16","PDrQ_bio17" ,"PWaQ_bio18",
#"PCoQ_bio19" ,"MDiR_bio2","Iso_bio3","TSe_bio4" ,"MTWaM_bio5",
#"MiTCoQ_bio6","TAR_bio7", "MTWeQ_bio8","MTDrQ_bio9","tmin","tmax","srad","vapr",
#"wind","elev","phh2o_15.30cm","nitrogen_15.30cm","soc_15.30cm","cec_15.30cm",
#"clay_15.30cm", "sand_15.30cm","silt_15.30cm","ocd_15.30cm",
#"ocs_0.30cm")

colnames(epafro_env1)= c("genotype","VINO","taxon","species", "lon","lat","country",
                         "AMT_1","MTWaQ_10","MTCoQ_11", "AP_12","PWeM_13","PDrM_14" ,"PSe_15",
                         "PWeQ_16","PDrQ_17" ,"PWaQ_18","PCoQ_19" ,"MDiR_2","Iso_3","TSe_4" ,
                         "MTWaM_5","MiTCoQ_6","TAR_7", "MTWeQ_8","MTDrQ_9","tmin","tmax","srad",
                         "vapr","wind","elev","phh2o","nitrogen","soc","cec","clay", "sand","silt",
                         "ocd","ocs")   

epafro_env1
str(epafro_env1)
write.csv(epafro_env1, 'epafro_allenvars.csv')
#epafro_env2 = epafro_env1
#head(epafro_env2)
#epafro_env2[epafro_env2$phh2o != "NA",]
epafro_wcl = epafro_env1[,c(1:32)]
epafro_soil = epafro_env1[,c(1:7,33:41)] 



library(factoextra)
library(FactoMineR)
epafro_pca.env=PCA(epafro_env1[,-c(1:8)], graph = F)
fviz_pca_biplot(epafro_pca.env)
epafro_env_cor = cor(epafro_env1[,-c(1:8)])
corrplot::corrplot(
  epafro_env_cor, order = "original", type = "upper", diag=T, tl.cex = 1,
  tl.col = c(rep("red", 23), "forestgreen", rep("blue", 9)), addCoef.col = "white",
  number.cex = 0.6, number.font = 1
)

#Calculation of the VIF for each environmental variable in the data set and 
# iterative pruning of the most collinear variables up to obtain a data set 
# with max(VIF) < vif.thr (i.e., the selected threshold for the VIF pruning)

# env_pruning()
# env_pruning applies a given VIF threshold to an environmental data set and 
# removes collinear variables in an iterative way up to obtain a pruned data set 
# with max(VIF) lower than the selected threshold.  
# 1. X: an environmental data set with individuals by row and the variables by
#       column
# 2. vif.thr: the selected VIF threshold
# Returns list with two slots:
# 1. the pruned environmental data set
# 2. the correlation matrix associated to the pruned data set 
env_pruning <- function(X = NULL, vif.thr = NULL) {
  
  cat("\n"); cat("When using this function, please cite:"); cat("\n")
  cat("[1] https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12372")
  cat("\n"); cat("\n")
  
  res <- list()
  
  env_vif <- multicol(vars = X)
  while (length(which(env_vif$VIF > vif.thr)) >= 1) {
    t <- rownames(env_vif)[1]
    t <- which(colnames(X) == t)
    X <- X[, -t]
    env_vif <- multicol(vars = X)
    print(env_vif); cat("\n")
  }
  
  res[[1]] <- X
  res[[2]] <- cor(X)
  names(res) <- c(paste0("Pruned environmental data set (VIF<", vif.thr, ")"), "Correlation matrix")
  return(res)
  
}



epafroenv_wcl_qced <- env_pruning(X = epafro_env1[,c(8:32)], vif.thr = 5)
epafroenv_soil_qced <- env_pruning(X = epafro_env1[,c(33:41)], vif.thr = 5)

# Now, we can compute the correlation matrix again and see if there is still 
# some 'pathological' correlation
corrplot(
  epafroenv_wcl_qced$`Correlation matrix`, 
  order = "original", type = "upper", diag=T, tl.cex = 1.2,
  tl.col = "black", 
  addCoef.col = "white",
  number.cex = 0.8, number.font = 1)

corrplot(
  epafroenv_soil_qced$`Correlation matrix`, 
  order = "original", type = "upper", diag=T, tl.cex = 1.2,
  tl.col = "black", 
  addCoef.col = "white",
  number.cex = 0.8, number.font = 1)



# Let's have a look to the distribution of the retained variables... 
hist_env(epafroenv_wcl_qced$`Pruned environmental data set (VIF<5)`)
hist_env(epafroenv_soil_qced$`Pruned environmental data set (VIF<5)`)

# Save the pruned environmental data set for future usage
write.table(epafroenv_wcl_qced$`Pruned environmental data set (VIF<5)`, 
            "epafroenv_wcl_qced.txt", col.names = TRUE, row.names = FALSE,
            sep = "\t", quote = FALSE)

write.table(epafroenv_soil_qced$`Pruned environmental data set (VIF<5)`, 
            "epafroenv_soil_qced.txt", col.names = TRUE, row.names = FALSE,
            sep = "\t", quote = FALSE)

#Final environmental variables
epafroenv_info_qced = epafro_env1[,c(1:7,13,17,18,25,29,33,34,37,39,40)]
head(epafroenv_info_qced)

write.csv(epafroenv_info_qced, 'epafro_envqced.csv', row.names = FALSE)
#allqced <- env_pruning(X = epafroenv_info_qced[,-c(1:7)], vif.thr = 5)

epafro_pca.env_qc=PCA(epafro_env1[,c(13,17,18,25,29,33,34,37,39,40)], graph = F)
fviz_pca_biplot(epafro_pca.env_qc)

library(plotly)
library(ggfortify)
library(ggplot2)

#PCA by species
epafrodf = epafro_env1[,c(13,17,18,25,29,33,34,37,39,40)]
pca_epafro_envall = prcomp(epafrodf, scale. = TRUE)
epafro_PCA_plot = autoplot(pca_epafro_envall, data = epafro_env1,  colour = "species", loadings = TRUE, loadings.colour = 'blue',
                           loadings.label = TRUE, loadings.label.size = 3) 
ggplotly(epafro_PCA_plot)

#PCA by country
epafrodf = epafro_env1[,c(13,17,18,25,29,33,34,37,39,40)]
pca_epafro_envall = prcomp(epafrodf, scale. = TRUE)
epafro_PCA_plot = autoplot(pca_epafro_envall, data = epafro_env1,  colour = "country", loadings = TRUE, loadings.colour = 'blue',
                           loadings.label = TRUE, loadings.label.size = 3) 
ggplotly(epafro_PCA_plot)

