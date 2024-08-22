##########################################
### da usare con environment : conda_py_39
##########################################
#
# This script easily help you resample/reproject your input data in a suitable format for FLEXTH
#
# The minimum input for FLEXTH is:
# 1. a binary (i.e. 1 - 0) raster map delineating flooded areas
# 2. topographical data in the form of a DTM (or a DEM if DTM is not available)
#
# Additional input may include:
# 3. a binary water body mask which delineates all permanent/seasonal water bodies
# 4. a binary no-data mask which delineates all areas where flood mapping was not performed
#
#
#
# The input raster files of FLEXTH must:
# 1. be in the same projected coordinate reference system
# 2. share the same extent and have the same resolution (i.e. the pixels must match)
#
#
# If you have a flood extent map in raster format (in a projected reference system) and a generic DTM for your area,
# this script resamples and reprojects your DTM to the extent, resolution and reference system of your input flood map. 
# Your flood map and the new DTM can then be used as inputs to FLEXTH.
#
# The script can be used with any other raster (not just with DTM). For example, it can be used for your water bodies and/or no-data masks.
# Just be careful to choose the appropriate resampling method by setting the variable "continuous_input" to True in case you are processing
# a continuous variable (DTM), or to False if you are processing a categorical variable (e.g. no-data or water bodies mask).
#
#
# INPUT:
# 1. input_flood_delineation:  you reference flood extent 
# 2. input_raster           :  the raster you want to resample/reproject 
#
# OUTPUT
# 1. output_raster          :  your processed input raster
#

#two versions: version 1 uses gdal and rasterio --- version 2 uses  only  gdal


#   INPUT FLOOD MAP RASTER MUST BE IN A PROJECTED REFERENCE SYSTEM


##########
#VERSION 1
##########



from osgeo import gdal
import glob, os
import logging
import rasterio
from pathlib import Path

gdal.UseExceptions()


# If you input file is a continuous variable (e.g. DTM), select "True".
# If it is a discrete (i.e. categorical) variable, select "False"

continuous_input = False


#SET INPUT-OUTUT
input_raster  = Path('/home/bettand/Desktop/storage_bettand/Documents/GFM_files_transfer/temp/congo/elaborazioni_mie/gfm/input/obswater.tif')
output_raster = Path('/home/bettand/Desktop/storage_bettand/Documents/GFM_files_transfer/temp/congo/elaborazioni_mie/gfm/input/obswater_processed_ZSTD.tif')


input_flood_delineation = Path('/home/bettand/Desktop/storage_bettand/Documents/GFM_files_transfer/temp/congo/elaborazioni_mie/gfm/input/flood.tif')




############
#DOES  STUFF
############


flood_dataset = rasterio.open(input_flood_delineation)
flood_array   = flood_dataset.read(1)

raster_dataset = rasterio.open(input_raster)


transform = flood_dataset.transform
flood_array_nrows, flood_array_ncol =  flood_array.shape
Dx= transform[0]
Dy= transform[4]
minX= transform[2]
maxY= transform[5]
maxX= minX + Dx*flood_array_ncol
minY= maxY + Dy*flood_array_nrows


input_proj       = raster_dataset.crs.wkt
output_proj      = flood_dataset.crs.wkt


if continuous_input == True:
    resampling_method = 'bilinear'# other options: 'average', 'cubic', 'lanczos' ... 
elif continuous_input == False:
    resampling_method = 'mode'   # other options: 'nearest' ...

gdal.Warp(str(output_raster), 
          str(input_raster), 
          outputBounds=[minX, minY, maxX, maxY],
          cropToCutline    = True,
          outputBoundsSRS = output_proj, 
          warpMemoryLimit= 5000,
          srcSRS= input_proj, 
          dstSRS= output_proj, 
          xRes=Dx,
          yRes=Dy, 
          resampleAlg= resampling_method,
          targetAlignedPixels = False,
          creationOptions = ["COMPRESS=ZSTD", "BIGTIFF=YES", "TILED=YES"]    # compression options: ZSTD, DEFLATE, LZW
          )




flood_dataset.close()
raster_dataset.close()





#############################
#######################################
#VERSION 2
######################################
#########################################################


from osgeo import gdal, osr, ogr
import os
from pathlib import Path
gdal.UseExceptions()


# If you input file is a continuous variable (e.g. DTM), select "True".
# If it is a discrete (i.e. categorical) variable, select "False"

continuous_input = True


#SET INPUT-OUTUT
input_raster  = Path('/home/bettand/Documents/GFM_files_transfer/temp/congo/elaborazioni_mie/congo_dtm_clip.tif')
output_raster = Path('/home/bettand/Documents/GFM_files_transfer/temp/congo/elaborazioni_mie/congo_reproject.tif')


input_flood_delineation = Path('/home/bettand/Documents/GFM_files_transfer/temp/congo/elaborazioni_mie/congo.tif')




flood_dataset_GDAL = gdal.Open(str(input_flood_delineation))
dem_dataset_GDAL = gdal.Open(str(input_raster))




geoTransform = flood_dataset_GDAL.GetGeoTransform()
minX = geoTransform[0]
maxY = geoTransform[3]
maxX = minX + geoTransform[1] * flood_dataset_GDAL.RasterXSize
minY = maxY + geoTransform[5] * flood_dataset_GDAL.RasterYSize

Dx= geoTransform[1]
Dy= geoTransform[5]

input_proj = dem_dataset_GDAL.GetSpatialRef()
output_proj = flood_dataset_GDAL.GetSpatialRef()




if continuous_input == True:
    resampling_method = 'bilinear'# other options: 'average', 'cubic', 'lanczos' ... 
elif continuous_input == False:
    resampling_method = 'mode'   # other options: 'nearest' ...




gdal.Warp(str(output_raster), 
          str(input_raster), 
          outputBounds=[minX, minY, maxX, maxY],
          outputBoundsSRS = output_proj, 
          warpMemoryLimit= 5000,
          srcSRS= input_proj, 
          dstSRS= output_proj, 
          xRes=Dx,
          yRes=-Dy, 
          resampleAlg= resampling_method,
          #targetAlignedPixels = True,
          creationOptions = ["COMPRESS=LZW", "BIGTIFF=YES", "TILED=YES"])


flood_dataset_GDAL = None
dem_dataset_GDAL   = None 
