devtools::install_github(repo = "alikhuseynov/seurat", ref = "feat/vizgen")
library(Seurat)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(magrittr)
library(BiocParallel)
library(progressr)
library(sf)

setwd("/kriegsteinlab/data3/LiWang/analysis/MERFISH/")

#load data
#ARKFrozen62PFC
ARKFrozen62PFC <- LoadVizgen(data.dir = "/kriegsteinlab/data3/LiWang/analysis/MERFISH/202208291531_20220829-Frozen62-GW24-section3_VMSC01801/region_0",  
                          fov = "ARKFrozen62PFC", 
                          assay = "Vizgen",
                          metadata = c("volume", "fov"), # add cell volume info
                          type = c("centroids"), # type of cell spatial coord matrices
                          z = 3L,
                          add.zIndex = TRUE, # add z slice section to a cell
                          update.object = TRUE,
                          use.BiocParallel = TRUE,
                          workers.MulticoreParam = 20, # for `BiocParallel` processing
                          #min.area = 5, # minimal polygon area to use as a threshold for filtering segmentaion geometries
                          add.molecules = TRUE, # if to add "molecules" coordinates to FOV of the object
                          verbose = TRUE
)

saveRDS(ARKFrozen62PFC, "analysis/ARKFrozen62PFC.rds")


#NIH4365BA10
NIH4365BA10 <- LoadVizgen(data.dir = "/kriegsteinlab/data3/LiWang/analysis/MERFISH/202208181215_20220818-4365-GW34_VMSC01801/region_0",  
                             fov = "NIH4365BA10", 
                             assay = "Vizgen",
                             metadata = c("volume", "fov"), # add cell volume info
                             type = c("centroids"), # type of cell spatial coord matrices
                             z = 3L,
                             add.zIndex = TRUE, # add z slice section to a cell
                             update.object = TRUE,
                             use.BiocParallel = TRUE,
                             workers.MulticoreParam = 20, # for `BiocParallel` processing
                             #min.area = 5, # minimal polygon area to use as a threshold for filtering segmentaion geometries
                             add.molecules = TRUE, # if to add "molecules" coordinates to FOV of the object
                             verbose = TRUE
)

saveRDS(NIH4365BA10, "analysis/NIH4365BA10.rds")



#UCSF2018003MFG
UCSF2018003MFG <- LoadVizgen(data.dir = "/kriegsteinlab/data3/LiWang/analysis/MERFISH/202209261319_20220926UCSF2018003repeat_VMSC01801/region_0",  
                            fov = "UCSF2018003MFG", 
                            assay = "Vizgen",
                            metadata = c("volume", "fov"), # add cell volume info
                            type = c("centroids"), # type of cell spatial coord matrices
                            z = 3L,
                            add.zIndex = TRUE, # add z slice section to a cell
                            update.object = TRUE,
                            use.BiocParallel = TRUE,
                            workers.MulticoreParam = 20, # for `BiocParallel` processing
                            #in.area = 5, # minimal polygon area to use as a threshold for filtering segmentaion geometries
                            add.molecules = TRUE, # if to add "molecules" coordinates to FOV of the object
                            verbose = TRUE
)

saveRDS(UCSF2018003MFG, "analysis/UCSF2018003MFG.rds")


ARKFrozen65V1 <- LoadVizgen(data.dir = "/kriegsteinlab/data3/LiWang/analysis/MERFISH/202307131325_20230713ARKFrozen65GW24V1_VMSC08202/region_0",  
                          fov = "ARKFrozen65V1", 
                          assay = "Vizgen",
                          metadata = c("volume", "fov"), # add cell volume info
                          type = c("centroids"), # type of cell spatial coord matrices
                          z = 3L,
                          add.zIndex = TRUE, # add z slice section to a cell
                          update.object = TRUE,
                          use.BiocParallel = TRUE,
                          workers.MulticoreParam = 12, # for `BiocParallel` processing
                          #min.area = 5, # minimal polygon area to use as a threshold for filtering segmentaion geometries
                          add.molecules = TRUE, # if to add "molecules" coordinates to FOV of the object
                          verbose = TRUE
)

saveRDS(ARKFrozen65V1, "analysis/ARKFrozen65V1.rds")

#NIH5900BA17
NIH5900BA17 <- LoadVizgen(data.dir = "/kriegsteinlab/data3/LiWang/analysis/MERFISH/202307201703_NIH-5900-BA17-GW34_VMSC08202/region_0",  
                          fov = "NIH5900BA17", 
                          assay = "Vizgen",
                          metadata = c("volume", "fov"), # add cell volume info
                          type = c("centroids"), # type of cell spatial coord matrices
                          z = 3L,
                          add.zIndex = TRUE, # add z slice section to a cell
                          update.object = TRUE,
                          use.BiocParallel = TRUE,
                          workers.MulticoreParam = 20, # for `BiocParallel` processing
                          #min.area = 5, # minimal polygon area to use as a threshold for filtering segmentaion geometries
                          add.molecules = TRUE, # if to add "molecules" coordinates to FOV of the object
                          verbose = TRUE
)

saveRDS(NIH5900BA17, "analysis/NIH5900BA17.rds")

#NIH4392BA17
NIH4392BA17 <- LoadVizgen(data.dir = "/kriegsteinlab/data3/LiWang/analysis/MERFISH/202310021545_20231002NIH-4392-BA17_VMSC13702/region_0",  
                          fov = "NIH4392BA17", 
                          assay = "Vizgen",
                          metadata = c("volume", "fov"), # add cell volume info
                          type = c("centroids"), # type of cell spatial coord matrices
                          z = 3L,
                          add.zIndex = TRUE, # add z slice section to a cell
                          update.object = TRUE,
                          use.BiocParallel = TRUE,
                          workers.MulticoreParam = 20, # for `BiocParallel` processing
                          #min.area = 5, # minimal polygon area to use as a threshold for filtering segmentaion geometries
                          add.molecules = TRUE, # if to add "molecules" coordinates to FOV of the object
                          verbose = TRUE
)

saveRDS(NIH4392BA17, "analysis/NIH4392BA17.rds")

