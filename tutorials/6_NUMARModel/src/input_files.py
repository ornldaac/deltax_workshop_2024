#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 10:07:12 2023

@author: alchrist
"""


import numpy as np
import rasterio
import os
import fnmatch
# import earthaccess
import urllib.request
from pathlib import Path


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


def get_input_files(AOI,aoi_dir,tmp_dir,EPSG,bounds_4326):

    Path(tmp_dir).mkdir(parents=True, exist_ok=True)

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
    ##########################################################################################
    ## AVIRIS-NG Aboveground Biomass maps
    ## Units are Mg/Ha
    print('')
    variable = 'agb'
    source = 'AVIRISNG'
    print('##### %s %s \n\n' %(source, variable))
    final_file = aoi_dir / ('%s_%s.tif' %(source,variable))

    if os.path.isfile(final_file)==False:
        print('# Search EarthData for %s %s over the AOI\n' %(source, variable))
        search = earthaccess.search_data(short_name = 'DeltaX_L3_AVIRIS-NG_AGB_V2_2138',
                                                bounding_box = tuple(bounds_4326),
                                                granule_name = '*%s*.tif' %(variable))
        print('# Download %s %s products: \n' %(source, variable))
        for i in range(len(search)):
            download_files = earthaccess.download(search[i], aoi_dir) 

        files = [os.path.join(dirpath,f)
            for dirpath,dirnames, files in os.walk(aoi_dir)
            for f in fnmatch.filter(files,'*%s*.tif' %(variable))]
        print('# Merging %s tiles\n' %(variable)) #EPSG:4326++5733
        with open("%s/%s.txt" %(tmp_dir,variable), 'w') as f:
            for item in files:
                f.write("%s\n" % item) 
        os.system("gdalbuildvrt %s/%s.vrt -separate -input_file_list %s/%s.txt -q"
                  %(tmp_dir,variable,tmp_dir,variable))
        vrt = rasterio.open('%s/%s.vrt' %(tmp_dir,variable))
        p = vrt.profile
        vrt = vrt.read()
        vrt = np.where(vrt==-9999,np.nan,vrt)
        avg = np.nanmean(vrt,axis=0)
        p['count']=1
        p['driver']='GTiff'
        with rasterio.open(tmp_dir/('%s_mean.tif' %(variable)),'w',**p) as dst:
            dst.write_band(1,avg.astype(float))

        print('# Crop the %s file to AOI and resample to  %sm \n' %(variable,res))
        os.system('gdalwarp -overwrite -t_srs epsg:%s -co COMPRESS=DEFLATE  %s/%s_mean.tif %s/%s.tif -q '%(EPSG,tmp_dir,variable,tmp_dir,variable))
        os.system('gdalwarp -overwrite -tr %s %s %s/%s.tif %s '\
            ' -te %s %s %s %s -srcnodata -9999 -dstnodata -9999 -co COMPRESS=DEFLATE -q'
            %(res,res,tmp_dir,variable,final_file,ulx,lry,lrx,uly))

        print('--> Saved as %s\n' %(final_file))
        print('')
    src = rasterio.open(final_file)
    AGB = np.where(src.read_masks(1)==0,np.nan,src.read(1))

    ##########################################################################################
    ##########################################################################################




    ##########################################################################################
    ##########################################################################################
    print('')
    variable = 'basins'
    source = 'HUC'
    print('##### %s %s \n\n' %(source, variable))
    final_file = aoi_dir / ('%s_%s.tif' %(source,variable))
    if os.path.isfile(final_file)==False:

        basin_file = aoi_dir / 'Delta-X_Basins_modified.tif'
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
    A = np.where(BASINS == 1,0,np.where(BASINS==2,1,np.nan))

    ##########################################################################################
    ##########################################################################################


    ##########################################################################################
    ##########################################################################################
    ## mineral deposition results from hydro/morphodynamic models
    ## units cm/year?
    print('')
    variable = 'IMAR'
    source = 'Delft3D'
    print('##### %s %s \n\n' %(source, variable))
    final_file = aoi_dir / ('%s_%s.tif' %(source,variable))
    if os.path.isfile(final_file)==False:


        print('# Search EarthData for Delft3D Atchafalaya Sediment Model Data over the AOI\n')
        search = earthaccess.search_data(short_name = 'DeltaX_Delft3D_Atchafalaya_MRD_2302',
                                    bounding_box = tuple(bounds_4326),
                                    granule_name = '*%s*' %(variable))
        print('# Download IMAR Delft3D data products: \n')
        download_files = earthaccess.download(search[0], aoi_dir) 

        files = [os.path.join(dirpath,f)
            for dirpath,dirnames, files in os.walk(aoi_dir)
            for f in fnmatch.filter(files,'*Atchafalaya*%s*.nc4' %(variable))][0]
        filename = files.split('/')[-1].split('.')[0]
        print('# Convert the Delft3D IMAR netCDF4 file to Geotiff\n')
        os.system('gdalwarp -q %s %s/%s.tif' %(files,tmp_dir,filename))
        os.system('gdalwarp -q -t_srs epsg:%s %s/%s.tif %s/%s_%s.tif' %(EPSG,tmp_dir, filename,  tmp_dir, filename,EPSG))
        
        print('# Search EarthData for Delft3D Terrebonne Sediment Model Data over the AOI\n')
        search = earthaccess.search_data(short_name = 'DeltaX_Delft3D_Terrebonne_MRD_2301',
                                    bounding_box = tuple(bounds_4326),
                                    granule_name = '*%s*' %(variable))
        print('# Download %s %s data products: \n' %(source, variable))
        download_files = earthaccess.download(search[0], aoi_dir) 

        files = [os.path.join(dirpath,f)
            for dirpath,dirnames, files in os.walk(aoi_dir)
            for f in fnmatch.filter(files,'*Terrebonne*%s*.nc4' %(variable))][0]
        filename = files.split('/')[-1].split('.')[0]
        print('# Convert the %s %s file to Geotiff\n' %(source, variable))
        os.system('gdalwarp -q %s %s/%s.tif' %(files,tmp_dir,filename))
        os.system('gdalwarp -q -t_srs epsg:%s %s/%s.tif %s/%s_%s.tif' %(EPSG,tmp_dir, filename,  tmp_dir, filename,EPSG))

        files = [os.path.join(dirpath,f)
            for dirpath,dirnames, files in os.walk(tmp_dir)
            for f in fnmatch.filter(files,'*%s*_%s.tif' %(variable,EPSG))]
        print('# Merging %s tiles\n' %(variable)) #EPSG:4326++5733
        with open("%s/%s.txt" %(tmp_dir,variable), 'w') as f:
            for item in files:
                f.write("%s\n" % item) 
        os.system("gdalbuildvrt %s/%s.vrt -separate -input_file_list %s/%s.txt -q"
                  %(tmp_dir,variable,tmp_dir,variable))
        vrt = rasterio.open('%s/%s.vrt' %(tmp_dir,variable))
        p = vrt.profile
        vrt = vrt.read()
        vrt = np.where(vrt==-9999,np.nan,vrt)
        avg = np.nanmean(vrt,axis=0)
        p['count']=1
        p['driver']='GTiff'
        with rasterio.open(tmp_dir/ ('%s_mean.tif'%(variable)),'w',**p) as dst:
            dst.write_band(1,avg.astype(float))

        print('# Crop the Delft3D IMAR file to AOI and resample to %sm \n' %(res))
        os.system('gdalwarp -overwrite -tr %s %s %s/%s_mean.tif %s '\
                ' -te %s %s %s %s -srcnodata -9999 -dstnodata -9999 -co COMPRESS=DEFLATE -q'
                %(res,res,tmp_dir,variable,final_file,ulx,lry,lrx,uly))
        print('--> Saved as %s\n' %(final_file))
        print('')


    IMAR = rasterio.open(final_file).read(1)
    ##########################################################################################
    ##########################################################################################

    ##########################################################################################
    ##########################################################################################
    ## new water mask
    ## 1 = water
    ## 0 = not water
    print('')
    variable = 'watermask'
    source = 'DeltaX'
    print('##### %s %s \n\n' %(source, variable))
    final_file = aoi_dir / ('%s_%s.tif' %(source,variable))

    if os.path.isfile(final_file)==False:
        print('# Search EarthData for %s %s over the AOI\n' %(source, variable))
        search = earthaccess.search_data(short_name = 'DeltaX_DEM_MRD_LA_2181',
                                                bounding_box = tuple(bounds_4326),
                                                granule_name = '*%s*.tif' %(variable))
        print('# Download %s %s products: \n' %(source, variable))
        for i in range(len(search)):
            download_files = earthaccess.download(search[i], aoi_dir) 

        files = [os.path.join(dirpath,f)
            for dirpath,dirnames, files in os.walk(aoi_dir)
            for f in fnmatch.filter(files,'*%s*.tif' %(variable))][0]
        print('# Crop the %s file to AOI and resample to  %sm \n' %(variable,res))
        os.system('gdalwarp -overwrite -t_srs epsg:%s -co COMPRESS=DEFLATE  %s %s/%s.tif -q '%(EPSG,files,tmp_dir,variable))
        os.system('gdalwarp -overwrite -tr %s %s %s/%s.tif %s '\
            ' -te %s %s %s %s -srcnodata -9999 -dstnodata -9999 -co COMPRESS=DEFLATE -q'
            %(res,res,tmp_dir,variable,final_file,ulx,lry,lrx,uly))

        print('--> Saved as %s\n' %(final_file))
        print('')
    WATER = rasterio.open(final_file).read(1)
    watermask = np.where(WATER ==0,1,np.nan)

    ##########################################################################################
    ##########################################################################################
    


    #########################################################################################
    #########################################################################################
    # Vegetation type from CPRA, which includes salinity
    # 1 = other/urban
    # 2 = Swamp
    # 3 = Fresh
    # 4 = Intermediate
    # 5 = Brackish
    # 6 = Saline
    # 7 = Permanent Water
    print('')
    variable = 'salinity'
    source = 'CPRA'
    print('##### %s %s \n\n' %(source, variable))
    final_file = aoi_dir / ('%s_%s.tif' %(source,variable))
    if os.path.isfile(final_file)==False:
        salinity_file = aoi_dir / 'vegtype2021'/'raster'/'2021_Veg_Zones_in_NLCD_CCAP_Mode_Raster_Data_031722.tif'
        if os.path.isfile(salinity_file)==False:
            if os.path.isfile(aoi_dir / 'vegtype2021.zip')==False:
                print('# Downloading CRPA Veg Zip file\n')
                urllib.request.urlretrieve('https://cims.coastal.louisiana.gov/Viewer/metadata/zips/vegtype2021.zip',aoi_dir / 'vegtype2021.zip')
            os.system('unzip -qq %s -d %s' %(aoi_dir / 'vegtype2021.zip',aoi_dir))
            os.remove(aoi_dir/'vegtype2021.zip')
        print('# Reprojecting\n')
        filename = '2021_Veg_Zones_in_NLCD_CCAP_Mode_Raster_Data_031722'
        os.system('gdalwarp -t_srs epsg:%s %s %s/%s_%s.tif -q' %(EPSG,salinity_file,tmp_dir,filename,EPSG))
        print('# Crop the Salinity file to AOI and resample to  %sm \n' %(res))
        os.system('gdalwarp -overwrite -tr %s %s %s/%s_%s.tif %s '\
                ' -te %s %s %s %s -srcnodata -9999 -dstnodata -9999 -co COMPRESS=DEFLATE -q'
                %(res,res,tmp_dir,filename,EPSG,final_file,ulx,lry,lrx,uly))
        print('--> Saved as %s\n' %(final_file))
        print('')

    VEGCLASS = rasterio.open(final_file).read(1)
    B = np.where( VEGCLASS == 3, 0,np.where((VEGCLASS==4) | (VEGCLASS ==5) , 1, np.where(VEGCLASS==6,2,np.nan)))

    ##########################################################################################
    ##########################################################################################
    CLASSES = ((A*3)+B)
    mask = np.where((watermask == 1) & (CLASSES <=6),1,np.nan)


    return AGB, VEGCLASS, IMAR, BASINS, CLASSES, mask, lats, lons


    
    












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


