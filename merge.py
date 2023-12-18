###################################################################################
## Creates the max flood extent by merging multiple rasters with flood delineations
###################################################################################



from rasterio.merge import merge
from rasterio.enums import Resampling 
import os


#input folder where the separater raster files are <<<<<<<<<
input_dir        =  '***'


#set output folder and file name <<<<<<<<
output_dir       =  '***'
output_file_name =  '***'




#list files in the input directory
input_paths     = os.listdir(input_dir)   
    
#full paths of files
input_list_full = [ os.path.join(input_dir,path)   for path in input_paths]
    

#mosaic rasters - select resampling method:  resampling=Resampling.**** with **** = {nearest, mode}
merge(input_list_full, nodata=255, resampling=Resampling.nearest, method='max', target_aligned_pixels=False, dst_path=os.path.join(output_dir,output_file_name))
    
    

    
    
    
    
    
####################################################################
## Resample and reprojects input raster in a given coordinate system
####################################################################



from rasterio.merge import merge
from rasterio.enums import Resampling 
import os
from osgeo import gdal


###set input-outputs <<<<<<<<
input_raster   =  '***'
output_raster  =  '***'


# select a PROJECTED reference system  <<<<<<<
output_proj      = 'ESRI:54009'   


# select resolution  <<<<<<<
Dx = 20 
Dy = 20


#warp! 
gdal.Warp(output_raster, 
          input_raster, 
          warpMemoryLimit= 5000,
          dstSRS= output_proj, 
          xRes=Dx,
          yRes=Dy, 
          resampleAlg='mode',
          targetAlignedPixels = False,
          creationOptions = ["COMPRESS=LZW", "BIGTIFF=YES", "TILED=YES"]
          )

