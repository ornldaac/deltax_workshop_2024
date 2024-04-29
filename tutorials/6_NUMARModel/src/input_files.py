#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 10:07:12 2023

@author: alchrist
"""


import numpy as np
import geopandas as gpd
import rasterio
from pathlib import Path
import os
import matplotlib.pyplot as plt
import xarray as xr
from shapely import Polygon
from tqdm.auto import tqdm  # provides a progressbar
import requests
import fnmatch
import earthaccess
import urllib.request


def calculate_belowground_biomass(AGB,mask, CLASSES,veg_e,bgb_agb):
    #b0  = 0.1154 * np.ones((height,width))        ## g cm-3	Bulk density of organic matter 
    BGB_AGB_ratio = bgb_agb[CLASSES]*mask
    e = veg_e[CLASSES]*mask
    r0 = ((AGB/1000*BGB_AGB_ratio)*e)/(1-np.exp(-e*50))
    # BGB = (AGB / 0.2) * landmask #(mg/ha)
    print('#########################')
    print('[[r0 ==> BELOWGROUND BIOMASS]]')
    # print('BGB = AGB * 0.2')
    print('BGB:AGB = 0.35 fresh marsh, 0.62 brackish marsh, and 1.26 saline marsh')
    print('BGB at surface (r0) = (r50 * e)/(1-exp(-e*50))')
    print('Attenuation rate (e) = 0.03 fresh marsh, 0.05 brackish marsh, and 0.0406 saline marsh')
    print('Units = g/cm2')
    print('#########################\n')    
    return r0*mask

def calculate_organicmatter_loading(si, CLASSES,mask,oms_si):
    oms = oms_si[CLASSES]*si*mask
    print('#########################')
    print('[[oms ==> ORGANIC MASS ACCUMULATION RATE]]')
    print('Units = g/cm2/yr')
    print('#########################\n')    

    return oms*mask


def get_input_files(AOI,aoi_dir,ref_dir,tmp_dir,EPSG,bounds_4326):
	
    # earthaccess.login()
    ## Input Files
    ## Should be EPSG 4326
    profile = AOI.profile
    transform = profile['transform']
    width = profile['width']
    height = profile['height']

    cols, rows = np.meshgrid(np.arange(width), np.arange(height))
    xs, ys = rasterio.transform.xy(AOI.transform, rows, cols)
    lats,lons = np.array(ys),np.array(xs)

    ulx = transform[2]
    uly = transform[5]
    lrx = ulx + width*transform[0]
    lry = uly + height*transform[4]
    res = profile['transform'][0]

    ##########################################################################################
    ## AVIRIS-NG Aboveground Biomass maps
    ## Units are Mg/Ha
    print('')
    print('##### ABOVEGROUND BIOMASS\n\n')
    final_file = aoi_dir / 'AVIRIS_AGB.tif'
    if os.path.isfile(final_file)==False:
        print('# Search EarthData for AVIRIS-NG Aboveground Biomass Data over the AOI\n')
        AGB_search = earthaccess.search_data(short_name = 'DeltaX_L3_AVIRIS-NG_AGB_V2_2138',
                                                bounding_box = tuple(bounds_4326),
                                                granule_name = '*.tif')
        print('# Download ABOVEGROUND BIOMASS from Spring 2021 and Fall 2021 AVIRIS-NG data products: \n')
        for i in range(len(AGB_search)):
            download_files = earthaccess.download(AGB_search[i], ref_dir) 

        agb_tifs = [os.path.join(dirpath,f)
            for dirpath,dirnames, files in os.walk(ref_dir)
            for f in fnmatch.filter(files,'*agb*.tif')]
        print('# Merging AGB tiles\n') #EPSG:4326++5733
        with open("%s/agb.txt" %(tmp_dir), 'w') as f:
            for item in agb_tifs:
                f.write("%s\n" % item) 
        os.system("gdalbuildvrt %s/agb.vrt -separate -input_file_list %s/agb.txt -q"
                  %(tmp_dir,tmp_dir))
        vrt = rasterio.open('%s/agb.vrt' %(tmp_dir))
        p = vrt.profile
        vrt = vrt.read()
        vrt = np.where(vrt==-9999,np.nan,vrt)
        avg = np.nanmean(vrt,axis=0)
        p['count']=1
        p['driver']='GTiff'
        with rasterio.open(tmp_dir/'agb_mean.tif','w',**p) as dst:
            dst.write_band(1,avg.astype(float))

        print('# Crop the WorldCover file to AOI and resample to  %sm \n' %(res))
        os.system('gdalwarp -overwrite -t_srs epsg:%s -co COMPRESS=DEFLATE  %s/agb_mean.tif %s/agb.tif -q '%(EPSG,tmp_dir,tmp_dir))
        os.system('gdalwarp -overwrite -tr %s %s %s/agb.tif %s '\
            ' -te %s %s %s %s -srcnodata -9999 -dstnodata -9999 -co COMPRESS=DEFLATE -q'
            %(res,res,tmp_dir,final_file,ulx,lry,lrx,uly))

        print('--> Saved as %s\n' %(final_file))
        print('')
    src = rasterio.open(final_file)
    AGB = np.where(src.read_masks(1)==0,np.nan,src.read(1))


    
    ##########################################################################################
 
    ##########################################################################################
    ## Landcover from World Cover (replace with AVIRIS-NG landcover when available)
    ## 10 = forest
    ## 20 = shrubland
    ## 30 = grassland
    ## 40 = cropland
    ## 50 = urban
    ## 60 = bare/sparse vegetation
    ## 70 = snow and ice
    ## 80 = permanent water body
    ## 90 = herbaceous wetland
    ## 95 = mangroves
    ## 100 = moss/lichen
    print('')
    print('##### ESA WorldCover landcover type\n\n')
    final_file = aoi_dir / 'ESA_WorldCover2021.tif'
    if os.path.isfile(final_file)==False:

        print('# Download ESA WorldCover 2021\n')    
        source = 'WorldCover'
        bounds = (bounds_4326[0]-1, bounds_4326[1]-1, bounds_4326[2]+1, bounds_4326[3]+1)
        geometry = Polygon.from_bounds(*bounds)

        ##Use AWS cloud data
        s3_url_prefix = "https://esa-worldcover.s3.eu-central-1.amazonaws.com"
        # load worldcover grid
        url = f'{s3_url_prefix}/v100/2020/esa_worldcover_2020_grid.geojson'
        grid = gpd.read_file(url)
        # get grid tiles intersecting AOI
        tiles = grid[grid.intersects(geometry)]
        for tile in tqdm(tiles.ll_tile):
            print('\t', tile)
            out_fn = ref_dir / f"ESA_WorldCover_10m_2021_v200_{tile}_Map.tif"
            if os.path.isfile(out_fn)==False:
                url = f"{s3_url_prefix}/v200/2021/map/ESA_WorldCover_10m_2021_v200_{tile}_Map.tif"
                r = requests.get(url, allow_redirects=True)
                with open(out_fn, 'wb') as f:
                    f.write(r.content)
        ### Authenticate to the Terrascope platform (registration required)
        ### create catalogue object and authenticate interactively with a browser
        # catalogue = Catalogue().authenticate()
        # products = catalogue.get_products("urn:eop:VITO:ESA_WorldCover_10m_2020_V1", geometry=geometry)
        # catalogue.download_products(products, folders[1],force = True)
        landcovers = [os.path.join(dirpath,f)
            for dirpath,dirnames, files in os.walk(ref_dir)
            for f in fnmatch.filter(files,'*_Map.tif')]
        with open("%s/%s.txt" %(tmp_dir,source), 'w') as f:
            for item in landcovers:
                f.write("%s\n" % item)
        # Merge NASADEM tiles to make topography file
        print('# Merging landcover tiles\n') #EPSG:4326++5733
        os.system("gdalbuildvrt %s/%s.vrt -input_file_list %s/%s.txt -a_srs EPSG:4326 -q"
                  %(tmp_dir,source,tmp_dir,source))
        print('# Crop the WorldCover file to AOI and resample to  %sm \n' %(res))
        os.system('gdalwarp -overwrite -t_srs epsg:%s -co COMPRESS=DEFLATE  %s/%s.vrt %s/%s.tif -q '%(EPSG,tmp_dir,source,tmp_dir,source))
        os.system('gdalwarp -overwrite -tr %s %s %s/%s.tif %s '\
            ' -te %s %s %s %s -srcnodata -9999 -dstnodata -9999 -co COMPRESS=DEFLATE -q'
            %(res,res,tmp_dir,source,final_file,ulx,lry,lrx,uly))
        print('--> Saved as %s\n' %(final_file))
        print('')

    LANDCOVER = rasterio.open(final_file).read(1)

    ##########################################################################################

    ##########################################################################################
    ## Vegetation type from CPRA, which includes salinity
    ## 1 = other/urban
    ## 2 = Swamp
    ## 3 = Fresh
    ## 4 = Intermediate
    ## 5 = Brackish
    ## 6 = Saline
    ## 7 = Permanent Water
    print('')
    print('##### CPRA Vegetation Type (with Salinity Zones)\n')
    final_file = aoi_dir / 'CPRA_Salinity.tif'
    if os.path.isfile(final_file)==False:

       
        salinity_file = ref_dir / 'vegtype2021'/'raster'/'2021_Veg_Zones_in_NLCD_CCAP_Mode_Raster_Data_031722.tif'
        if os.path.isfile(salinity_file)==False:
            if os.path.isfile(ref_dir / 'vegtype2021.zip')==False:
                print('# Downloading CRPA Veg Zip file\n')
                urllib.request.urlretrieve('https://cims.coastal.louisiana.gov/Viewer/metadata/zips/vegtype2021.zip',ref_dir / 'vegtype2021.zip')
            os.system('unzip -qq %s -d %s' %(ref_dir / 'vegtype2021.zip',ref_dir))
            os.remove(ref_dir/'vegtype2021.zip')
        print('# Reprojecting\n')
        filename = '2021_Veg_Zones_in_NLCD_CCAP_Mode_Raster_Data_031722'
        os.system('gdalwarp -t_srs epsg:%s %s %s/%s_%s.tif -q' %(EPSG,salinity_file,tmp_dir,filename,EPSG))
        print('# Crop the Salinity file to AOI and resample to  %sm \n' %(res))
        os.system('gdalwarp -overwrite -tr %s %s %s/%s_%s.tif %s '\
                ' -te %s %s %s %s -srcnodata -9999 -dstnodata -9999 -co COMPRESS=DEFLATE -q'
                %(res,res,tmp_dir,filename,EPSG,final_file,ulx,lry,lrx,uly))
        print('--> Saved as %s\n' %(final_file))
        print('')

    SALINITY = rasterio.open(final_file).read(1)

    ##########################################################################################
    



    ##########################################################################################
    print('')
    print('##### Modified Hydrologic Basins\n')
    final_file = aoi_dir / 'HUC_basins.tif'
    if os.path.isfile(final_file)==False:

        basin_file = ref_dir / 'Delta-X_Basins_modified.tif'
        print('# Reprojecting\n')
        filename = 'Delta-X_Basins_modified'
        os.system('gdalwarp -t_srs epsg:%s %s %s/%s_%s.tif -q' %(EPSG,basin_file,tmp_dir,filename,EPSG))
        print('# Crop the Basin file to AOI and resample to  %sm \n' %(res))
        os.system('gdalwarp -overwrite -tr %s %s %s/%s_%s.tif %s '\
                ' -te %s %s %s %s -srcnodata -9999 -dstnodata -9999 -co COMPRESS=DEFLATE -q'
                %(res,res,tmp_dir,filename,EPSG,final_file,ulx,lry,lrx,uly))
        print('--> Saved as %s\n' %(final_file))
        print('')

    BASINS = rasterio.open(final_file).read(1)
   
    ##########################################################################################

    A = np.where(BASINS == 1,0,np.where(BASINS==2,1,np.nan))
    B = np.where( SALINITY == 3, 0,np.where((SALINITY==4) | (SALINITY ==5) , 1, np.where(SALINITY==6,2,np.nan)))
    CLASSES = ((A*3)+B)

    ##########################################################################################
    ## mineral deposition results from hydro/morphodynamic models
    ## units cm/year?
    print('')
    print('##### Delft3D Inorganice Mass Accumulation Rate (IMAR)\n\n')
    final_file = aoi_dir / 'Delft3D_IMAR.tif'
    if os.path.isfile(final_file)==False:


        print('# Search EarthData for Delft3D Atchafalaya Sediment Model Data over the AOI\n')
        IMAR_search = earthaccess.search_data(short_name = 'DeltaX_Delft3D_Atchafalaya_MRD_2302',
                                    bounding_box = tuple(bounds_4326),
                                    granule_name = '*IMAR*')
        print('# Download IMAR Delft3D data products: \n')
        download_files = earthaccess.download(IMAR_search[0], ref_dir) 

        morph_file = [os.path.join(dirpath,f)
            for dirpath,dirnames, files in os.walk(ref_dir)
            for f in fnmatch.filter(files,'*IMAR*.nc4')][0]
        filename = morph_file.split('/')[-1].split('.')[0]
        print('# Convert the Delft3D IMAR netCDF4 file to Geotiff\n')
        os.system('gdalwarp -q %s %s/%s.tif' %(morph_file,tmp_dir,filename))
        os.system('gdalwarp -q -t_srs epsg:%s %s/%s.tif %s/%s_%s.tif' %(EPSG,tmp_dir, filename,  tmp_dir, filename,EPSG))
        
        print('# Search EarthData for Delft3D Terrebonne Sediment Model Data over the AOI\n')
        IMAR_search = earthaccess.search_data(short_name = 'DeltaX_Delft3D_Terrebonne_MRD_2301',
                                    bounding_box = tuple(bounds_4326),
                                    granule_name = '*IMAR*')
        print('# Download IMAR Delft3D data products: \n')
        download_files = earthaccess.download(IMAR_search[0], ref_dir) 

        morph_file = [os.path.join(dirpath,f)
            for dirpath,dirnames, files in os.walk(ref_dir)
            for f in fnmatch.filter(files,'*IMAR*.nc4')][0]
        filename = morph_file.split('/')[-1].split('.')[0]
        print('# Convert the Delft3D IMAR netCDF4 file to Geotiff\n')
        os.system('gdalwarp -q %s %s/%s.tif' %(morph_file,tmp_dir,filename))
        os.system('gdalwarp -q -t_srs epsg:%s %s/%s.tif %s/%s_%s.tif' %(EPSG,tmp_dir, filename,  tmp_dir, filename,EPSG))

        imar_tifs = [os.path.join(dirpath,f)
            for dirpath,dirnames, files in os.walk(tmp_dir)
            for f in fnmatch.filter(files,'*IMAR*_%s.tif' %(EPSG))]
        print('# Merging IMAR tiles\n') #EPSG:4326++5733
        with open("%s/imar.txt" %(tmp_dir), 'w') as f:
            for item in imar_tifs:
                f.write("%s\n" % item) 
        os.system("gdalbuildvrt %s/imar.vrt -separate -input_file_list %s/imar.txt -q"
                  %(tmp_dir,tmp_dir))
        vrt = rasterio.open('%s/imar.vrt' %(tmp_dir))
        p = vrt.profile
        vrt = vrt.read()
        vrt = np.where(vrt==-9999,np.nan,vrt)
        avg = np.nanmean(vrt,axis=0)
        p['count']=1
        p['driver']='GTiff'
        with rasterio.open(tmp_dir/'imar_mean.tif','w',**p) as dst:
            dst.write_band(1,avg.astype(float))

        print('# Crop the Delft3D IMAR file to AOI and resample to %sm \n' %(res))
        os.system('gdalwarp -overwrite -tr %s %s %s/imar_mean.tif %s '\
                ' -te %s %s %s %s -srcnodata -9999 -dstnodata -9999 -co COMPRESS=DEFLATE -q'
                %(res,res,tmp_dir,final_file,ulx,lry,lrx,uly))
        print('--> Saved as %s\n' %(final_file))
        print('')


    IMAR = rasterio.open(final_file).read(1)
    ##########################################################################################

    ##########################################################################################
    ## new water mask
    ## 1 = water
    ## 0 = not water
    print('##### Delta-X Watermask')
    final_file = aoi_dir/'watermask.tif'
    if os.path.isfile(final_file)==False:

        watermask_file = ref_dir  / 'Delta-X_watermask.tif'
        os.system('gdalwarp -q -t_srs epsg:%s %s %s' %(EPSG,watermask_file,str(watermask_file)[:-4] + '_%s.tif' %(EPSG)))
        os.system('gdalwarp -overwrite -tr %s %s %s %s '\
                ' -te %s %s %s %s -srcnodata -9999 -dstnodata -9999 -co COMPRESS=DEFLATE -q'
                %(res,res,str(watermask_file)[:-4] + '_%s.tif' %(EPSG),final_file,ulx,lry,lrx,uly))
    WATER = rasterio.open(aoi_dir/ 'watermask.tif').read(1)
    watermask = np.where(WATER ==0,1,np.nan)
    ##########################################################################################
    
    mask = np.where((watermask == 1) & (CLASSES <=6),1,np.nan)


    return AGB, SALINITY, LANDCOVER, IMAR, BASINS, CLASSES, mask, lats, lons


    
    












# ## defaults
# b   = 2.5000 * np.ones((height,width))        ## g cm-2	Total aboveground biomass
# b0  = 0.1154 * np.ones((height,width))        ## g cm-3	Bulk density of organic matter 
# c0  = 0.1500 * np.ones((height,width))        ## g g-1	Lignin content in leaf litter
# c1  = 0.1000 * np.ones((height,width))        ## g g-1	Ash concentration in litter
# e   = 0.0400 * np.ones((height,width))        ## cm-1	A distribution parameter for root biomass
# f1  = 0.7000 * np.ones((height,width))        ## g g-1	Ratio of litter leaves to total litter
# f2  = 0.4500 * np.ones((height,width))        ## g g-1	Proportion of microbial respiration during decomposition
# f3  = 0.0004 * np.ones((height,width))        ## g g-1	Proportion of labile organic matter flowing into refractory organic matter after decomposition
# fc1	= 0.2500 * np.ones((height,width))        ## g g-1	Fraction of fine root that is refractory
# fc2	= 0.3000 * np.ones((height,width))        ## g g-1	Fraction of large root that is refractory
# ka  = 0.9000 * np.ones((height,width))        ## y-1	Decay constant for labile organic matter on the surface
# kb  = 0.2560 * np.ones((height,width))        ## y-1	Decay constant for labile organic matter below the suface
# kc  = 0.0010 * np.ones((height,width))        ## y-1	Decay constant for refractory organic matter
# ke  = 0.4000 * np.ones((height,width))        ## y-1	Litter export rate
# km  = 0.0425 * np.ones((height,width))        ## y-1	Turnover rate of wood
# kr  = 0.1000 * np.ones((height,width))        ## y-1	Turnover rate of fine roots
# kt  = 0.2760 * np.ones((height,width))        ## y-1	Decay constant of litter twigs
# kw  = 0.0830 * np.ones((height,width))        ## y-1	Decay constant of dead wood
# lp  = 0.1300 * np.ones((height,width))        ## g cm-2 y -1	Litter Production
# r0  = 0.0740 * np.ones((height,width))        ## g cm-2	Total root biomass at the surface
# si  = 0.1950 * np.ones((height,width))        ## g cm-2 y-1	Annual deposition of mineral sediment on the surface


