# FLEXTH

FLEXT is a tool to enhance flood maps (e.g. satellite-derived) by accounting for terrain topography. It expands inundations to areas which are likely to be flooded based on their altimetry and provides estimates of water levels and water depths. The algorithm requires, as a primary input, a flood delineation map and a DTM. Additional information may include areas excluded from flood mapping and/or permanent water bodies. All input must be provided via georeferenced raster format (GeoTIFF) in a suitable projected reference system. An additional script is provided named "DTM_2_floodmap.py" which can easily  help you resample and reproject your input DTM into the same grid, extent and reference system as your flood map raster. 
Since water depth is the primary proxy for flood damages, the tool aims to facilitate flood impact assessment over large scales with minimum supervision and quick computational times.



CITE :  A.Betterle & P.Salamon - Water depth estimate and flood extent enhancement for satellite-based inundation maps  - ***



For further information contact: 

andrea.betterle@ec.europa.eu  
peter.salamon@ec.europa.eu
