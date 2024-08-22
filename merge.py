###################################################################################
## This script can be used to merge a set of input binary raster files (i.e 0-1).
## For example, it can generate the maximum flood extent from a set of flood delineation files.
## It can also be used for other layers (e.g. no-data or permanent water dodies)
###################################################################################



from rasterio.merge import merge
from rasterio.enums import Resampling 
import rasterio
import os, glob
from pathlib import Path

os.environ['GDAL_CACHEMAX'] = '8192'



#input folder where the separater raster files are <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< INPUT
input_dir        =  Path('/home/bettand/Desktop/storage_bettand/Documents/GFM_files_transfer/temp/congo/exclusion_obswater/') 



#set output directory  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< INPUT
output_dir       =  Path('/home/bettand/Desktop/storage_bettand/Documents/GFM_files_transfer/temp/congo/exclusion_obswater/') 



# choose output name <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< INPUT
output_file_name =  'obswater.tif'
#output_file_name =  'exclusion.tif'


##########################
## CHOOSE ONE:  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< INPUT
    ##########


###
#1 : all files in a given directory 
###

#list files in the input directory
input_paths     = os.listdir(input_dir)   
  
#full paths of files
input_list_full = [ os.path.join(input_dir,path)   for path in input_paths]
    
####
## 2: all files in subdirectories having the "layer_to_merge" string in their name   
####

#layer_to_merge = 'EXCLAYER'
#layer_to_merge = 'REFERENCE_WATER_OUT'
layer_to_merge = 'OBSERVED_WATER'
#layer_to_merge = 'EXCLUSION_LAYER'
input_list_full = glob.glob( str(input_dir/f"**/*{layer_to_merge}*.tif")  , recursive=True)



##########
## MOSAICS
##########

#mosaic rasters - select resampling method:  resampling=Resampling.**** with **** = {nearest, mode}
merge(input_list_full, 
      dst_path=os.path.join(output_dir,output_file_name),
      nodata=255, 
      resampling=Resampling.nearest,
      method='max', 
      target_aligned_pixels=False
      )
    



#####################################################################
 ########################################################################
######################################  
    
    
    
####################################################################
## Resample and reprojects input raster in a given coordinate system
####################################################################


# if a reference raster "reference_raster" is given, uses that and matches the pixels
# otherwise explicitly defines target CRS and resolution



from rasterio.merge import merge
from rasterio.enums import Resampling 
import os
from osgeo import gdal
import rasterio
from pathlib import Path

gdal.UseExceptions()

###set input-outputs <<<<<<<<
input_raster      =  Path('/home/bettand/Documents/GFM_files_transfer/brazil/RioGrandeDolSol_28April-14May2024/exclusion/exclusion.tif')
reference_raster  =  Path('/home/bettand/Documents/GFM_files_transfer/brazil/RioGrandeDolSol_28April-14May2024/50m/input/flood.tif')

output_raster     =   Path('/home/bettand/Documents/GFM_files_transfer/brazil/RioGrandeDolSol_28April-14May2024/50m/input/exclusion.tif')



try:
    reference_raster
    print("Reference data provided")
    reference     =  rasterio.open( reference_raster )
    output_proj   =  str(reference.crs)
    trans         = reference.transform
    Dx = trans[0] 
    Dy = abs(trans[4])
    outputBounds  =  reference.bounds
except NameError:
    print("Reference data not provided")
    # select a PROJECTED reference system  <<<<<<<
    output_proj      = 'ESRI:54009'   
    outputBounds     = None
    # select resolution  <<<<<<<
    Dx = 20 
    Dy = 20





#warp! 
gdal.Warp(output_raster, 
          input_raster, 
          outputBounds = outputBounds,
          warpMemoryLimit= 5000,
          dstSRS= output_proj, 
          xRes=Dx,
          yRes=Dy, 
          resampleAlg='mode',
          targetAlignedPixels = False,
          creationOptions = ["COMPRESS=LZW", "BIGTIFF=YES", "TILED=YES"]
          )

