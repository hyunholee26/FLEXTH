# Welcome to FLEXTH - the flood extent enhancement and water depth estimation 
# tool for (satellite-derived) inundation maps

# The tool enhances flood delineation maps (e.g. satellite-derived) by extending 
# floods based on terrain topography and additional optional information. 
# Furthermore, the algorithm employs topographical data (in the form of a DTM) 
# in combination with flood delineations to provide water depth estimates across 
# the flooded regions. The algorithm requires, as a primary input, a flood delineation
# map and a DTM. Additional information may include areas excluded from flood mapping 
# and/or permanent water bodies. All input must be provided via binary georeferenced raster 
# format (GeoTIFF) in a suitable projected reference system.
# Since water depth is the primary proxy for flood damages, the tool aims to facilitate
# flood impact assessment over large scales with minimum supervision and quick
# computational times.
# 
# The script 'DTM_2_floodmap.py' helps to easily bring your DTM into the same extent, 
# resolution, grid and projections  as the input flood map raster.
#
# For further details see Betterle and Salamon (2024) on Natural Hazards and 
# Earth System Sciences, or contact andrea.betterle@ec.europa.eu or peter.salamon@ec.europa.eu



import numpy as np
import rasterio
import cv2
from astropy.convolution import convolve, Tophat2DKernel, Box2DKernel
import warnings
from scipy.spatial.distance import cdist
from scipy import signal
import os
import time
from pathlib import Path
from scipy.spatial import cKDTree
from skimage import morphology, measure




################
###################
######################
#INPUTs and PARAMETERS 
#########################
############################
###############################



#   input_dir must contain the following maps in a suitable PROJECTED reference system :
#
#   flood.tif        :   binary map delineating flooded areas
#   dtm.tif          :   DTM
#   exclusion.tif    : * binary map delineating areas excluded from flood mapping
#   permanent_water  : * binary map delineating permanent/seasonal water bodies
#   obswater         : * binary map delineating all observed water (union between flood and perm_water) - to be used if permanent_water.tif is not provided
#   (* denotes optional input)
#
#  output
#  WD.tif : expanded flood extent with water depth estimates (in cm)
#  WL.tif : expanded flood extent with water level estimates (in m a.s.l.)
#
# All variable with names starting with "param" denote adjustable parameters.
# The most relevant are listed below whereas throughout the codes there are 
# additional less sensitive and more operative ones. Testing showed that the 
# default values are effective and robust in a wide range of settings. 
# Nonetheless, parameters can be tweaked to match specific needs and/or use cases.



# set working directories
input_dir  = Path( '/.../..../input' )
output_dir = Path( '/.../..../output' )

 
#Water level estimation method (options: 'method_A', 'method_B'): 
param_WL_estimation_method = 'method_A'  


# Select output: Water depth ("WD") , Water level ("WL"), both ("WL_WD")
param_output_map = "WL_WD"

#parameters
param_threshold_slope        =  0.05   # S_max : pixels steeper than this are not used to estimate water level
param_size_gaps_close        =  0.05   # A_g   : up to this size (in km2) gaps in the flood map will be closed

param_border_percentile      =  0.50    # P: assign water level based on the percentile P of the elevation of border cells (valid if Method B is selected)
param_max_number_neighbors   =  100     # N_max: number of border pixels used to compute water levels for a single continguos flooded area
param_min_flood_size         =  100     # N_min: if the number of valid border pixels along the border is less than this, WLs are estimated based on the distribution of pixels inside the flooded area
param_inner_quantile         =  0.98    # P*: for  flooded areas that don't meet the criteria above, uses this percentile of the elevation of flooded pixels to estimate water levels
param_inverse_dist_exp       =  2       # alpha: inverse distance weighting exponent used to interpolate WL inside flooded areas 

#param_threshold_distance_factor_1  =  0.5    # a: factor that regulates flood expansion
#param_threshold_distance_factor_2  =  0.5    # b: exponent that regulates flood expansion

param_max_propagation_distance     =  10     # maximum propagation distance in km
param_distance_range               =  100    # flooded areas of this size (km2) reach half of the maximum propagation distance  

param_WD_star  =  10  # WD* dummy water depth in cm assigned when estimated WD in flooded areas are negative



####################
###########################################################################
## no need to modify anything after this point (unless you have good ideas)
###########################################################################
#####################################



############
##############
#FUNCTIONS######
##################
####################

    
#  identifies the indices of neighboring cells in an array 
#  i,j : indices of the target location; 
#  n_row and n_col:  number or rows and columns of the  target array 
#  connectivity: type of connectivity to use (4 , 8)
 
def ij_neighbors(i,j,n_row,n_col,connectivity):  
    if connectivity == 8:
        #8d-connected
        neighbors=np.array([[i-1,j-1],[i-1,j], [i-1,j+1], [i,j+1], [i+1,j+1], [i+1,j], [i+1,j-1], [i,j-1]])  #8-connected
        for indice in range(len(neighbors)):
            if neighbors[indice,0]> n_row-1 or neighbors[indice,0]<0 :
                neighbors[indice,0]=i
                neighbors[indice,1]=j
            if neighbors[indice,1]> n_col-1 or neighbors[indice,1]<0 :
                neighbors[indice,0]=i
                neighbors[indice,1]=j    
         
    elif connectivity == 4:
        #4d-connected
        neighbors=np.array([ [i-1,j], [i,j+1], [i+1,j],  [i,j-1]])  
        for indice in range(len(neighbors)):
            if neighbors[indice,0]> n_row-1 or neighbors[indice,0]<0 :
                neighbors[indice,0]=i
                neighbors[indice,1]=j
            if neighbors[indice,1]> n_col-1 or neighbors[indice,1]<0 :
                neighbors[indice,0]=i
                neighbors[indice,1]=j                     
                
    return neighbors.astype('uint16')     
   
    


# function that computes the weighted quantiles out of a PDF

def weighted_quantile(values, percentile, weights=None):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: 
    :computes one percentile across the columns of a 2d array
    :percentile should be in [0, 1]!
    :param values: numpy.array with data
    :param percentile: the value of the quantile
    :param sample_weight: array-like of the same length as `array`
    :return: numpy.array with computed percentile for each row of data.
    """
    values = np.array(values)
    percentile = np.array(percentile)
    if weights is None:
        weights= np.ones(np.shape(values))
    weights = np.array(weights)
    assert np.all(percentile >= 0) and np.all(percentile <= 1), \
        'percentile should be in [0, 1]'


    sorter = np.argsort(values)
    values = np.take_along_axis(values, sorter, axis = 1)   
    weights = np.take_along_axis(weights, sorter, axis = 1)      

    weighted_quantiles = np.cumsum(weights, axis = 1 ) - 0.5 * weights

    weighted_quantiles /= np.sum(weights, axis = 1)[:, np.newaxis] * np.ones ( np.shape(weighted_quantiles))
    
    
    return values [np.array( np.arange(0,np.shape(values)[0],1)  ), np.argmin( (weighted_quantiles - percentile  )**2, axis = 1) ] 




# =============================================================================
# # Move masked values to the right in each row
# def justify_2d(arr, mask_value):
#     mask  =  arr!=mask_value
#     shifted = np.array([np.concatenate((row[mask[i]], row[~mask[i]])) for i, row in enumerate(arr)])
#     return shifted
# 
# 
# =============================================================================

###############
# main function 
 ##############

def flood_processing():
       
    flood_georeferenced  =   rasterio.open(  input_dir / 'flood.tif' )
    flood                =   rasterio.open(  input_dir / 'flood.tif' ).read(1)
    dem                  =   rasterio.open(  input_dir / 'dtm.tif'   ).read(1)
    
    
    
    
    #   checks if optional inputs are provided
    # - exclusion mask : areas masked from flood mapping
    # - observed water : flood water + permanent waters
    # - permanent water: permanent and/or seasonal water bodies
    
    if os.path.isfile(input_dir / 'exclusion.tif'):
        exclusion      =  rasterio.open(  input_dir / 'exclusion.tif' ).read(1)
    else:
        exclusion      = np.full_like(flood, 0)
        
        
    if os.path.isfile(input_dir / 'obswater.tif'):
        obswater      =  rasterio.open(  input_dir / 'obswater.tif' ).read(1)
    else:
        obswater      =  np.copy(flood)
    
    
    if os.path.isfile(input_dir / 'permanent_water.tif'):
        permanent_water      =  rasterio.open(  input_dir / 'permanent_water.tif' ).read(1)
    else:
        permanent_water      =  np.full_like(flood, 0)
    
    
    
    
    #geotranformation
    trans =  flood_georeferenced.transform
    
    trans*(0,0)  # goes from COL-ROW index to X-Y:  trans*(ncol,nrow)=(x_coordinate,y_coordinate)
    
    #extract the dimension of the raster
    size_equi7_tile= dem.shape    
    n_row, n_col = size_equi7_tile 
    
    # pixel size in m
    L = np.abs(np.array(trans*(0,1))[1]-np.array(trans*(0,0))[1])  
     
    #local slope
    slope_max =  np.gradient( dem, L)
    slope_max =  np.sqrt(slope_max[0]**2 + slope_max[1]**2 )   
    

    slope_mild_binary   =   np.array(slope_max <  param_threshold_slope).astype('uint8')
    
    
    
    
    #############################################################################
    ##connected components anaysis with cv2 - identifies contiguous flooded areas        
    #############################################################################
    
    #makes sure there are no novalues
    exclusion[exclusion == 255]               = 0  
    flood[flood == 255]                       = 0
    obswater[obswater == 255]                 = 0
    permanent_water[permanent_water == 255]   = 0
    
     
    
    #choose kernel 
    kernel_1       = np.ones((3,3),np.uint8)
    kernel_1_cross = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype = np.uint8)
    kernel_2 = np.ones((5,5),np.uint8)
    
    
    #crates the binary map to be used as input for the connected components analysis
    flood_binary                  = np.copy(flood).astype('uint8')
    flood_binary[flood_binary>0]  = 1  
    
    #runs morphological closing to remove small holes and very irregular borders
    flood_binary    = cv2.morphologyEx(flood_binary    ,    cv2.MORPH_CLOSE, kernel_1_cross,  iterations = 2 )
    
    
    
    
    #same as above
    exclusion_binary  = np.copy(exclusion).astype('uint8')
    exclusion_binary[exclusion_binary>0] =1
    exclusion_binary = cv2.dilate(exclusion_binary  ,  kernel_1)
    
    
    #same as above
    obswater_binary   = np.copy(obswater).astype('uint8')
    obswater_binary[obswater_binary>0] = 1
    obswater_binary = cv2.morphologyEx(obswater_binary ,    cv2.MORPH_CLOSE, kernel_1_cross,  iterations = 2 )
    
    
    ################################
    #closes small holes in flood map
    ################################
    
    print ('closing small gaps in flood map...')
    threshold_n_pixels = (param_size_gaps_close*1e6 / L**2).astype(int)
    
    binary_map_complementary_flood    = np.logical_not(flood_binary)

    # Label connected components and calculate their areas
    labeled_map, num_features = measure.label(binary_map_complementary_flood, return_num=True)
    component_areas = measure.regionprops(labeled_map)

    # Iterate over the component areas and perform morphological operations based on the threshold
    for component in component_areas:
        if component.area < threshold_n_pixels:
            flood_binary    [labeled_map == component.label]    = 1  # Fill small gaps
            obswater_binary [labeled_map == component.label] = 1




    
    
    
    
    
    if os.path.isfile(input_dir / 'permanent_water.tif'):
        #same as above
        permanent_water_binary                    = np.copy(permanent_water).astype('uint8')
        permanent_water_binary[permanent_water_binary>0] = 1
        permanent_water_binary = cv2.morphologyEx(permanent_water_binary ,    cv2.MORPH_CLOSE, kernel_1_cross,  iterations = 2 )
        
    else:
        permanent_water_binary = obswater_binary - flood_binary
        permanent_water_binary[(permanent_water_binary!=0) & (permanent_water_binary!=1) ] = 0
    

    
    del exclusion, obswater, permanent_water, binary_map_complementary_flood
    
    
    
    # !!!  opencv with connectivity = 4 has problems in some version of the package !!!!
           
    # Applying cv2.connectedComponents()
    # possibly change the connectivity type 
    # z_stats: [xtop, ytop, xwidth, ywidth, #pixels] - first element of stats correspondes to background value
    
    z_num_labels, z_labels, z_stats, z_centroids = cv2.connectedComponentsWithStats(flood_binary,connectivity=8)   #<---------  
    
      
    flood_binary_dilated          =   cv2.dilate(flood_binary  ,  kernel_1)
    flood_binary_eroded           =   cv2.erode(flood_binary   ,  kernel_1_cross)
    
    exclusion_binary_dilated      =   cv2.dilate(exclusion_binary ,  kernel_1)

    permanent_water_binary_dilated =   cv2.dilate(permanent_water_binary ,  kernel_1)
    
    
    #flood_dilated_label = np.copy(z_labels)
    flood_dilated_label           =   cv2.dilate(z_labels.astype('float32')  ,  kernel_1).astype('uint32')

    
    flood_border_binary           =   flood_binary_dilated    - flood_binary_eroded 
    flood_border_binary_inner     =   flood_binary            - flood_binary_eroded
    
    
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    flood_border_label            =   flood_dilated_label    
    #flood_border_label            =    z_labels     #- flood_eroded_label
    
    
    slope_mild_binary_dilated      =   cv2.dilate(slope_mild_binary  ,  kernel_1)
    
    
    exclusion_binary_dilated      [ (exclusion_binary       ==0)  & (exclusion_binary_dilated       ==1)  & (flood_binary==0) & (flood_binary_dilated==1)] = 0
    permanent_water_binary_dilated[ (permanent_water_binary ==0)  & (permanent_water_binary_dilated ==1)  & (flood_binary==0) & (flood_binary_dilated==1)] = 0
    slope_mild_binary_dilated     [ (slope_mild_binary      ==0)  & (slope_mild_binary_dilated      ==1)  & (flood_binary==0) & (flood_binary_dilated==1)] = 0
    
    
    flood_border_binary           =   (  flood_border_binary  *  abs(exclusion_binary_dilated.astype('int8') - 1)   *  abs(permanent_water_binary_dilated.astype('int8') - 1) ).astype('uint8')   *    slope_mild_binary_dilated  
    
    dem_border_flood              =   dem * flood_border_binary
      
    
    flood_border_binary           =  flood_border_binary * flood_border_binary_inner
          
    
    ## AVERAGES THE DEM ALONG THE BORDER 

    dem_border_flood_smooth= convolve ( dem_border_flood, kernel_1, mask = flood_border_binary == 0, preserve_nan = True )
    dem_border_flood_smooth[flood_border_binary==0]=np.nan     
    
    

    ######################################################################
    ##FLOOD DEPTH ESTIMATION IN INITIALLY FLOODED AREAS
    ######################################################################
    
    
    z_water_level = np.full(size_equi7_tile, np.NaN, dtype='float32')
    
    ROW, COL   =  np.mgrid[0:n_row,0:n_col].astype('uint32')
    
    row_flood                          =   np.ma.masked_array(ROW,      z_labels==0).compressed()
    col_flood                          =   np.ma.masked_array(COL,      z_labels==0).compressed()
    flood_labels_compressed            =   np.ma.masked_array(z_labels, z_labels==0).compressed()
        
    row_flood_border                   =   np.ma.masked_array(ROW,                     flood_border_binary==0).compressed()
    col_flood_border                   =   np.ma.masked_array(COL,                     flood_border_binary==0).compressed()        
    flood_borders_labels_compressed    =   np.ma.masked_array(flood_border_label,      flood_border_binary==0).compressed()
    dem_border_flood_smooth_compressed =   np.ma.masked_array(dem_border_flood_smooth, flood_border_binary==0).compressed()

    
    # assigns water levels to initial flood areas based on the elevation of the 
    # edges of flooded areas that don't border with the exclusion mask

     
    print('Assigns water elevation to initial flooded areas')
    

    for i in range (1,len(z_stats)):
        print (f'percent progress:  {int( i/len(z_stats) *100 ) }')
        

        temp_position_flood           =  np.stack((row_flood[flood_labels_compressed==i],col_flood[flood_labels_compressed==i]),axis = 1)
        temp_position_border          =  np.stack((row_flood_border[flood_borders_labels_compressed==i],col_flood_border[flood_borders_labels_compressed==i]),axis = 1)
        temp_dem_border_flood_smooth  =  dem_border_flood_smooth_compressed[flood_borders_labels_compressed==i]
        
        len_border  = len(temp_position_border)
 

        if  len_border < param_min_flood_size: 
                
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')
                z_water_level[temp_position_flood[:,0], temp_position_flood[:,1] ]  =  np.nanquantile(dem[temp_position_flood[:,0], temp_position_flood[:,1] ], param_inner_quantile)
    
    
        else:   
    
            
            # Build a cKDTree for each set of points
            #tree_FLOOD  = cKDTree(temp_position_flood)
            tree_BORDER = cKDTree(temp_position_border)
    
    
            distances, indices = tree_BORDER.query(temp_position_flood, k=param_max_number_neighbors,workers=10)


            distances= distances.astype('float32')
            distances[distances==0]  = 1
            indices  = indices.astype('uint32')


            if distances.ndim == 1:
                distances = np.expand_dims(distances, axis=1)
                indices   = np.expand_dims(indices, axis=1)

            # =============================================================================
            # # You can also find all points in set2 within a certain distance from a point in set1
            # # For example, find all points in set2 within a radius of 2 from a point in set1
            # neighbor_indices = tree_P.query_ball_point(points_A, r=10)
            # =============================================================================


            if param_max_number_neighbors > len_border:
                indices   = indices  [:, :len_border]
                distances = distances[:, :len_border]


            
            #####################
            #####################
            ## choose the weights
            #####################
            #####################
                
                
            ################################
            # inverse distance weighting IDW
            ################################
                
            weights  =   1 / (distances ** param_inverse_dist_exp).astype('float32')

                
              
                
            ####################################
            # exponential distance weighting EDW
            ####################################

            #weights  =  np.exp ( -distances.astype * L * np.log(2) / param_exponential_IDW   )
                
               
                
             ##############################
             ##############################
             ## choose WL estimation method
             ##############################
             ##############################


            if param_WL_estimation_method == 'method_A':
            #METHOD A : water level is interpoleated using IDW or EDW weights
                    
                
                dem_temp_statistic    =  (np.sum((temp_dem_border_flood_smooth[indices] * weights).astype('float32'), axis = 1 )  / np.sum( weights , axis = 1 ) ).astype('float32')
                
                
            elif param_WL_estimation_method == 'method_B':
                #METHOD B : water level is a distence-weighted QUANTILE of all closest dem neighboring cells 
        
                    
                dem_temp_statistic    =   weighted_quantile(temp_dem_border_flood_smooth[indices], param_border_percentile, weights  =  weights)
                             
                    
    
    
            #plotta
            if False:
                min_ =np.min(temp_dem_border_flood_smooth) 
                max_ = np.max(temp_dem_border_flood_smooth)
                import matplotlib.pyplot as plt
                plt.scatter(temp_position_flood[:,0], temp_position_flood[:,1], c=dem_temp_statistic, cmap='viridis', marker='s')
                plt.clim(min_, max_)
                plt.scatter(temp_position_border[:,0], temp_position_border[:,1], c=temp_dem_border_flood_smooth, cmap='viridis', marker='s')
                plt.clim(min_, max_)
    
                
    
    
                
            z_water_level[temp_position_flood[:,0], temp_position_flood[:,1]]    =    dem_temp_statistic 
    
    
     
    

    ###############
    #################
    ##FLOOD EXPANSION
    ###################
    #####################
 
    
    print('Augmenting flooded areas')
    

    param_connectivity                 =  8         # <<INPUT connectivity used to spread flooded areas

     
 
    z_water_level_masked_compressed         =   np.ma.masked_array(z_water_level, flood_border_binary_inner==0).compressed() 
    z_water_level_masked_compressed_initial =   np.ma.masked_array(z_water_level, flood_border_binary_inner==0).compressed()
    row_flood_border                        =   np.ma.masked_array(ROW,           flood_border_binary_inner==0).compressed()
    col_flood_border                        =   np.ma.masked_array(COL,           flood_border_binary_inner==0).compressed()
    flood_labels_border                     =   np.ma.masked_array(z_labels,      flood_border_binary_inner==0).compressed()
    
    FLOOD = np.array ( [flood_labels_border,  row_flood_border,  col_flood_border,  z_water_level_masked_compressed, z_water_level_masked_compressed_initial ]  ).T
    
    #sorts based on decreasing water elevation
    FLOOD_sort_descending= np.flipud(FLOOD[FLOOD[:, 3].argsort()])
    
    #size of initial flooded areas in m2
    dim_flooded_area = z_stats[:,4] * L**2
      
    
    
    ########################
    # max expansion distance  
    ########################
   
    ##POWER LAW
    #threshold_distance = param_threshold_distance_factor_1 * np.power(dim_flooded_area, param_threshold_distance_factor_2 )
    
    ##EXPONENTIAL
    threshold_distance = param_max_propagation_distance * 1000   *  ( 1 -   pow(2,   -  dim_flooded_area   / (param_distance_range * 1e6  )  )        )
    
    
    #avoids singularities in case no flood expansion is set
    threshold_distance = threshold_distance + 1
    
    
        
    FLOOD_sort_descending_list= list (FLOOD_sort_descending)
      
    z_water_level_augmented = np.copy(z_water_level)
    
    
    z_labels_augmented      = np.copy(z_labels)
    temp_dist_from_border   = np.zeros(z_labels.shape, dtype = 'float16')
    
    
    
    i=0
    i_counter  =  len ( FLOOD_sort_descending_list )
    i_stop     =  len ( FLOOD_sort_descending_list )
    total_processed = 0
    
    
    
       
    while i < i_counter: 
        
        if i == i_stop:
            FLOOD_sort_descending = np.vstack(FLOOD_sort_descending_list)
            FLOOD_sort_descending = FLOOD_sort_descending[i_stop: ,:]
            FLOOD_sort_descending = np.flipud(FLOOD_sort_descending[FLOOD_sort_descending[:, 3].argsort()]) 
            FLOOD_sort_descending_list= list (FLOOD_sort_descending)
            i=0
            i_stop    = len(FLOOD_sort_descending)
            i_counter = len(FLOOD_sort_descending) 
    
    
        vicini =  ij_neighbors(FLOOD_sort_descending[i,1], FLOOD_sort_descending[i,2] , n_row, n_col, param_connectivity)  
         
     
    
        label     = FLOOD_sort_descending[i,0].astype('uint32')
        
        
        initial_wl= FLOOD_sort_descending[i,4]  
        

        temp_= ( L +  temp_dist_from_border[FLOOD_sort_descending[i,1].astype('uint16'), FLOOD_sort_descending[i,2].astype('uint16')]  ) / threshold_distance[label]
        
        wl = (  initial_wl  -  ( initial_wl  - dem[vicini[:,0],vicini[:,1]]   )  * temp_  )     * np.heaviside(1-temp_ , 0 )   +    dem[vicini[:,0],vicini[:,1]]  *     (  1- np.heaviside(1-temp_ , 0 )    )
       
        wl[wl>FLOOD_sort_descending[i,3]] = FLOOD_sort_descending[i,3]
        
        
                                                                                                                # 1. must be nodata (i.e. not alreay assigned)
        condition =  (   
            
                     np.heaviside( (exclusion_binary[vicini[:,0],vicini[:,1]] > 0)  + (obswater_binary[vicini[:,0],vicini[:,1]] > 0) + (permanent_water_binary[vicini[:,0],vicini[:,1]] > 0)   , 0).astype('bool')          # 2. must be in the exlusion or permanent water layer

                   * ( dem[vicini[:,0],vicini[:,1]] <  (wl-0.01) ) # FLOOD_sort_descending[i,3] #  z_flooded_areas_augmented[i][2]                   # 3. dem must be lower than the corresponding water level    
                 
                   * np.heaviside(  np.isnan( z_water_level_augmented[vicini[:,0],vicini[:,1]])         , 0  ).astype('bool') 

                     )
    
       
        temp_dist_from_border[vicini[:,0][condition],vicini[:,1][condition]]    = L * np.array([1.414, 1, 1.414, 1, 1.414, 1, 1.414, 1])[condition]     + temp_dist_from_border[FLOOD_sort_descending[i,1].astype('uint16'), FLOOD_sort_descending[i,2].astype('uint16')]
        
        z_water_level_augmented[vicini[:,0][condition],vicini[:,1][condition]]  =  wl[condition]   
        
        FLOOD_sort_descending_list.append  (  np.array([ label*(condition[condition]) ,  vicini[:,0][condition],  vicini[:,1][condition]  ,  wl[condition] , initial_wl*(condition[condition])  ]).T  ) 
        
        z_labels_augmented[vicini[:,0][condition],vicini[:,1][condition]]  =  label 
        
        i_counter+= np.sum(condition)
                    
                             
        i+=1 
        
        #print(i)
        total_processed +=1
    

    
    
    ###############
    ##smoothing WL
    ###############

     
    print('Smoothing \n')
           
    
    ## first closes small holes in the flood map and then smooths (only the extended flooded areas)
    
    param_kernels_smoothing_size = 5  # 5 ; 10
    param_number_smoothings      = 20  # 20;  1 
    

    
    #ker  = np.ones((param_kernels_smoothing_size,param_kernels_smoothing_size)) / param_kernels_smoothing_size**2
    #kernel = np.outer(signal.windows.gaussian(70, 8), signal.windows.gaussian(70, 8))
    
    ker  = np.array(Tophat2DKernel( int(param_kernels_smoothing_size/2) )).astype('float32')
    ker  = ker/ np.sum(ker)
    
    
    
    z_water_level_augmented_initial = np.copy(z_water_level_augmented)
    


    
    param_kern_closing_size     = 3 
    param_closing_iterations    = 2  
    
    
    kercloseSizes = (param_kern_closing_size, param_kern_closing_size)
    kernelclose   = cv2.getStructuringElement(cv2.MORPH_RECT, kercloseSizes)  #kernel shapes : cv2.MORPH_ELLIPSE; cv2.MORPH_RECT
    big_holes   = cv2.morphologyEx((~np.isnan(z_water_level_augmented_initial)).astype('uint8'), cv2.MORPH_CLOSE, kernelclose ,iterations = param_closing_iterations)
    big_holes   =  ~(big_holes.astype('bool'))
    
    
    z_water_level_augmented[ np.isnan(z_water_level_augmented_initial) ] = dem[  np.isnan(z_water_level_augmented_initial)  ]
    
    for v in range(param_number_smoothings):

        z_water_level_augmented  =  signal.convolve2d(z_water_level_augmented, ker, mode='same', boundary='symm', fillvalue=0)
    
        z_water_level_augmented[ ~np.isnan(z_water_level) ] = z_water_level[  ~np.isnan(z_water_level)  ]
        
        z_water_level_augmented[ (np.isnan(z_water_level_augmented_initial)) & (big_holes == 1 ) ] = dem[  np.isnan(z_water_level_augmented_initial)  & (big_holes == 1 ) ]
        
    
    
    
    z_water_level_augmented[ np.isnan(z_water_level_augmented_initial) ] = np.nan
    
    
        

     
    ######
    # WD #
    ######
    WD                        = ( ( z_water_level_augmented - dem ) * 100  ).astype('float32')   
    
    WD[WD<0] = 0 
    WD[WD > 0]                = WD[WD > 0] + param_WD_star
    WD[(WD==0) & (flood_binary>0)  ] = param_WD_star
    
    
    
    #dummy water depth assigned to permanent water bodies
    WD[permanent_water_binary==1] = 999 
    
    WD[np.isnan(WD)] = 0
    
    
    
    ######
    # WL #
    ######
 
    WL                        = ( z_water_level_augmented   ).astype('float32')   
    
    WL[WD > 0]                = WL[WD > 0] + param_WD_star/100
    WL[(WD==0) & (flood_binary>0)  ] = dem[(WD==0) & (flood_binary>0) ]  + param_WD_star/100
    
    
    
    #dummy water level assigned to permanent water bodies
    WL[permanent_water_binary==1] = 999 
    
    
    
    
    ########
    ##SAVING
    ########
    
    print('Saving \n')
    
    if param_output_map == "WD" or  param_output_map == "WL_WD" :
    
        with rasterio.open(
            output_dir / f'WD_{param_WL_estimation_method}_Smax_{param_threshold_slope}_Nmax_{param_max_number_neighbors}_a_{param_inverse_dist_exp}_Dmax_{param_max_propagation_distance}_A12_{param_distance_range}.tif',
            mode="w",
            driver="GTiff",
            compress='lzw',
            height=n_row,
            width=n_col,
            count=1,
            dtype=rasterio.uint16,
            crs=flood_georeferenced.crs,
            transform=flood_georeferenced.transform,
            nodata= 0
            ) as water_depth:
                water_depth.write(WD, 1)
                
                
    if param_output_map == "WL" or  param_output_map == "WL_WD" :
    
        with rasterio.open(
            output_dir / f'WL_{param_WL_estimation_method}_Smax_{param_threshold_slope}_Nmax_{param_max_number_neighbors}_a_{param_inverse_dist_exp}_Dmax_{param_max_propagation_distance}_A12_{param_distance_range}.tif',
            mode="w",
            driver="GTiff",
            compress='lzw',
            height=n_row,
            width=n_col,
            count=1,
            dtype=rasterio.float32,
            crs=flood_georeferenced.crs,
            transform=flood_georeferenced.transform,
            ) as water_level:
                water_level.write(WL, 1)
    
   
    
    print('Finish \n')
      
      
    
    
    flood_georeferenced.close()
    

        
########        
## RUN ! 
########     
  
start = time.time()     
    
if __name__ == '__main__':
    
    flood_processing()
    
        

end = time.time()


print(f'Elapsed time:{int((end - start)/60)} minutes')
