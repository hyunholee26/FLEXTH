##########################################
### da usare con environment : conda_py_39
##########################################

# script that takes in input a raster flood map and a dtm. 
# it brings the dtm to the same extent and resolution (with matching pixels) as the flood map


#two versions: version 1 uses gdal and rasterio --- version 2 uses  only  gdal


#  !!! INPUT FLOOD MAP RASTER MUST BE IN A PROJECTED REFERENCE SYSTEM


##########
#VERSION 1
##########



from osgeo import gdal
import glob, os
import logging
import rasterio
from pathlib import Path


#SET INPUT-OUTUT
input_DTM  = Path('/***/***/dtm_in.tif')
output_DTM = Path('/***/***/dtm_out.tif')


input_flood_delineation = Path('/***/***/flood.tif')




###############################################################################
#DOES ITS STUFF

flood_dataset = rasterio.open(input_flood_delineation)
flood_array   = flood_dataset.read(1)

dtm_dataset = rasterio.open(input_DTM)


transform = flood_dataset.transform
flood_array_nrows, flood_array_ncol =  flood_array.shape
Dx= transform[0]
Dy= transform[4]
minX= transform[2]
maxY= transform[5]
maxX= minX + Dx*flood_array_ncol
minY= maxY + Dy*flood_array_nrows


input_proj       = dtm_dataset.crs.wkt
output_proj      = flood_dataset.crs.wkt


gdal.Warp(str(output_DTM), 
          str(input_DTM), 
          outputBounds=[minX, minY, maxX, maxY],
          cropToCutline    = True,
          outputBoundsSRS = output_proj, 
          warpMemoryLimit= 5000,
          srcSRS= input_proj, 
          dstSRS= output_proj, 
          xRes=Dx,
          yRes=Dy, 
          resampleAlg='bilinear',
          targetAlignedPixels = False,
          creationOptions = ["COMPRESS=LZW", "BIGTIFF=YES", "TILED=YES"]
          )




flood_dataset.close()
dtm_dataset.close()






##########
#VERSION 2
##########



from osgeo import gdal, osr, ogr
import os
from pathlib import Path

#SET INPUT-OUTUT
input_DTM  = Path('/***/***/dtm_in.tif')
output_DTM = Path('/***/***/dtm_out.tif')


input_flood_delineation = Path('/***/***/flood.tif')



flood_dataset_GDAL = gdal.Open(input_flood_delineation)
dem_dataset_GDAL = gdal.Open(input_DTM)




geoTransform = flood_dataset_GDAL.GetGeoTransform()
minX = geoTransform[0]
maxY = geoTransform[3]
maxX = minX + geoTransform[1] * flood_dataset_GDAL.RasterXSize
minY = maxY + geoTransform[5] * flood_dataset_GDAL.RasterYSize

Dx= geoTransform[1]
Dy= geoTransform[5]

input_proj = dem_dataset_GDAL.GetSpatialRef()
output_proj = flood_dataset_GDAL.GetSpatialRef()




gdal.Warp(str(output_DTM), 
          str(input_DTM), 
          outputBounds=[minX, minY, maxX, maxY],
          outputBoundsSRS = output_proj, 
          warpMemoryLimit= 5000,
          srcSRS= input_proj, 
          dstSRS= output_proj, 
          xRes=Dx,
          yRes=-Dy, 
          resampleAlg='bilinear',
          #targetAlignedPixels = True,
          creationOptions = ["COMPRESS=LZW", "BIGTIFF=YES", "TILED=YES"])


flood_dataset_GDAL = None
dem_dataset_GDAL   = None 
