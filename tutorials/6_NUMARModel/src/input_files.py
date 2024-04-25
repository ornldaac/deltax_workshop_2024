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


# veg_c0 = np.array([0.30,0.07,0.02,0.06,0.03,0.06])
# veg_c1 = np.array([0,0,0,0,0,0])
# veg_c2 = np.array([0,0,0,0,0,0])
# veg_c4 = np.array([0.25,0.25,0.25,0.25,0.25,0.25])
# veg_e  = np.array([0.0276, 0.05, 0.0406, 0.0043, 0.02, 0.0345])
# veg_fc1 = np.array([0.24, 0.48, 0.41, 0.15, 0.41, 0.43])
# veg_kb  = np.array([0.90, 0.37, 0.37, 0.74, 0.58, 0.31])
# veg_kc  = np.array([0.0109, 0.0109, 0.0109, 0.0109, 0.0109, 0.0109])
# veg_kl  = np.array([0.39,0.39,0.39,0.39,0.39,0.39])
# veg_kr  = np.array([0.04, 2.65, 2.02, 0.64, 0.94, 0.58])
# bgb_agb = np.array([0.35, 0.62, 1.26, 0.35, 0.62, 1.26])
# oms_si = np.array([0.03,0.05,0.0406,0.03,0.05,0.0406])
# veg_b0 = np.array([0.04,0.0502,0.0554,0.0850,0.0560,0.0850])
# veg_bi = np.array([1.99, 1.99, 1.99, 1.99, 1.99, 1.99])

## Parameters based on AVIRIS biomass
def calculate_aboveground_biomass(AGB,mask):
    #b   = 2.5000 * np.ones((height,width))        
    ## g cm-2	Total aboveground biomass
    AGB_mgHA = np.where(AGB==-9999,np.nan,AGB)
    AGB_mgHA = AGB_mgHA * mask#(mg/ha)
    print('#########################')
    print('[[b ==> ABOVEGROUND BIOMASS]]')
    print('Convert units (AVIRIS-NG units are Mg/Ha)')
    print('Units = g/cm2')
    print('#########################\n')

    AGB_gCM2 = AGB_mgHA * 1000000 / 100000000
    #plt.imshow(AGB_gCM2)
    
    return AGB_gCM2
def calculate_belowground_biomass(AGB,mask, CLASSES,veg_e,bgb_agb):
    #b0  = 0.1154 * np.ones((height,width))        ## g cm-3	Bulk density of organic matter 
    BGB_AGB_ratio = bgb_agb[CLASSES]
    e = veg_e[CLASSES]
    r0 = ((AGB*BGB_AGB_ratio)*e)/(1-np.exp(-e*50))
    # BGB = (AGB / 0.2) * landmask #(mg/ha)
    print('#########################')
    print('[[r0 ==> BELOWGROUND BIOMASS]]')
    # print('BGB = AGB * 0.2')
    print('BGB:AGB = 0.35 fresh marsh, 0.62 brackish marsh, and 1.26 saline marsh')
    print('BGB at surface (r0) = (r50 * e)/(1-exp(e*50))')
    print('Attenuation rate (e) = 0.03 fresh marsh, 0.05 brackish marsh, and 0.0406 saline marsh')
    print('Units = g/cm2')
    print('#########################\n')    
    return r0*mask

def calculate_organicmatter_loading(si, CLASSES,mask,oms_si):
    oms = oms_si[CLASSES]*si
    print('#########################')
    print('[[oms ==> ORGANIC MASS ACCUMULATION RATE]]')
    print('Units = g/cm2/yr')
    print('#########################\n')    

    return oms*mask


# def calculate_e_distribution(BGB,landmask):
#     e = 0.04 * BGB * landmask
#     print('#########################')
#     print('[[e ==> DISTRIBUTION PARAMETER FOR ROOT BIOMASS]]')
#     print('e = 0.04')
#     print('Units = cm-1')
#     print('#########################\n')    
#     return e

# ## Parameters based on landcover type
def calculate_bulkdensity_OM(CLASSES,mask):
    #BD = np.where(landcover ==1, 1, 0.1154)
    print('#########################')
    print('[[b0 ==> BULK DENSITY OF ORGANIC MATTER]]')
    # print('BD = 0.1154')
    print('This value is based on salinity and sediment zone based on data from Castañeda-Moya and Solohin, 2023')
    print('Units = g/cm3')
    print('#########################\n')    
    # BD = 0.1154 * landcover * landmask
    b0 = veg_b0[CLASSES]
    return b0*mask

def calculate_bulkdensity_IM(CLASSES,mask):
    #BD = np.where(landcover ==1, 1, 0.1154)
    print('#########################')
    print('[[bi ==> BULK DENSITY OF INORGANIC MATTER]]')
    # print('BD = 0.1154')
    print('This uniform value based on data from Castañeda-Moya and Solohin, 2023')
    print('Units = g/cm3')
    print('#########################\n')    
    # BD = 0.1154 * landcover * landmask
    bi = veg_bi[CLASSES]
    return bi*mask

def calculate_lignin_leaflitter(CLASSES,mask):
    #lignin = np.where(landcover == 1,1,0.15)
    print('#########################')
    print('[[c0 ==> LIGNIN CONTENT IN LEAF LITTER]]')
    print('c0 = 0.15')
    print('Units = g lignin/g leaf litter')
    print('#########################\n')    
    # lignin = 0.15 *landcover * landmask
    c0 = veg_c0[CLASSES]
    return c0*mask

def calculate_ash_root(CLASSES,mask):
    #ash = np.where(landcover==1,1,0.1)
    print('#########################')
    print('[[c1 ==> ASH CONCENTRATION IN LITTER]]')
    print('c1 = 0.1')
    print('Units = g ash/ g litter')
    print('#########################\n')   
    # ash = 0.1 * landcover * landmask
    c1 = veg_c1[CLASSES]
    return c1*mask
# def calculate_ratio_leaflitter(landcover,landmask):
#     #ratio = np.where(landcover == 1,1,0.7)
#     print('#########################')
#     print('[[f1 ==> RATIO OF LEAF LITTER TO TOTAL LITTER]]')
#     print('f1 = 0.7')
#     print('Units = g leaf litter / g total litter')
#     print('#########################\n')    
#     ratio = 0.7 * landcover * landmask
#     return ratio

# def calculate_refractory_fineroot(landcover,landmask):
#     print('#########################')
#     print('[[fc1 ==> FRACTION REFRACTORY OF FINE ROOT]]')
#     print('fc1 = 0.25')
#     print('Units = g refractory / g fine root')
#     print('#########################\n')    
#     fine_refrac = 0.25 * landcover * landmask
#     return fine_refrac

def calculate_lignin_roots(CLASSES,mask):
    print('#########################')
    print('[[fc2 ==> LIGNIN FRACTION OF BELOWGROUND BIOMASS]]')
    print('fc1 = 0.3')
    print('Units = g refractory / g large root')
    print('#########################\n')    
    # fc1 = 0.3 * landcover * landmask
    fc1 = veg_fc1[CLASSES]
    return fc1*mask

# def calculate_lignin_BGB(landcover,landmask):
#     print('#########################')
#     print('[[fc1 ==> LIGNIN FRACTION OF BELOWGROUND BIOMASS]]')
#     print('fc1 = 0.25')
#     print('Units = g refractory / g fine root')
#     print('#########################\n')    
#     fc0 = 0.25 * landcover * landmask
#     return fc0

# def calculate_refractory_largeroot(landcover,landmask):
#     print('#########################')
#     print('[[fc2 ==> FRACTION REFRACTORY OF LARGE ROOT ]]')
#     print('fc2 = 0.2')
#     print('Units = g refractory / g large root')
#     print('#########################\n')    
#     large_refrac = 0.2 * landcover * landmask
#     return large_refrac

## Paramers based on salinity
# def calculate_microbial_resp_decomp(input_dst,landmask):
#     print('#########################')
#     print('[[f2 ==> PROPORTION OF MICROBIAL RESPIRATION DURING DECOMPOSITION]]')
#     print('f2 = 0.45')
#     print('Units = g/g')
#     print('#########################\n')    
#     respiration = 0.45 * input_dst * landmask
#     return respiration

# def calculate_labile2refractory_decomp(input_dst,landmask):
#     print('#########################')
#     print('[[f3 ==> PROPORTION OF LABILE OM FLOWING INTO REFRACTORY OM AFTER DECOMP ]]')
#     print('f3 = 0.004')
#     print('Units = g/g')
#     print('#########################\n')        
#     flow = 0.004 * input_dst * landmask
#     return flow

# def calculate_decay_labileOM_surface(input_dst,landmask):
#     print('#########################')
#     print('[[ka ==> DECAY CONSTANT FOR LABILE OM ON SURFACE ]]')
#     print('ka = 0.9')
#     print('Units = yr-1')
#     print('#########################\n')    
#     decay = 0.9 * input_dst * landmask
#     return decay

def calculate_decay_labileOM_below(CLASSES,mask):
    print('#########################')
    print('[[kb ==> DECAY CONSTANT FOR LABILE OM BELOW SURFACE ]]')
    print('kb = 0.256')
    print('Units = yr-1')
    print('#########################\n')        
    # kb = 0.256 * input_dst * landmask
    kb = veg_kb[CLASSES]
    return kb*mask

def calculate_decay_celluloseOM(CLASSES,mask):
    print('#########################')
    print('[[kc ==> DECAY CONSTANT FOR REFRACTORY OM ]]')
    print('kc = 0.001')
    print('Units = yr-1')
    print('#########################\n')    
    # kc = 0.001 * input_dst * landmask
    kc = veg_kc[CLASSES]
    return kc*mask

def calculate_litter_decomp(CLASSES,mask):
    print('#########################')
    print('[[kl ==> LIGNIN DECOMPOSITION RATE ]]')
    print('kl = 0.001')
    print('Units = yr-1')
    print('#########################\n')   
    kl = veg_kl[CLASSES]
    return kl*mask

# def calculate_litter_export(input_dst,landmask):
#     print('#########################')
#     print('[[ke ==> LITTER EXPORT RATE]]')
#     print('ke = 0.4')
#     print('Units = yr-1')
#     print('#########################\n')    
#     ke = 0.4 * input_dst * landmask
#     return ke

# def calculate_turnover_wood(input_dst,landmask):
#     print('#########################')
#     print('[[km ==> TURNOVER RATE OF WOOD]]')
#     print('km = 0.0425')
#     print('Units = yr-1')
#     print('#########################\n')    
#     km = 0.0425 * input_dst * landmask
#     return km

def calculate_fineroot_turnover(CLASSES,mask):
    print('#########################')
    print('[[kr ==> TURNOVER RATE OF FINE ROOTS]]')
    print('kr = 0.1')
    print('Units = yr-1')
    print('#########################\n')    
    # kr = 0.1 * input_dst * landmask
    kr = veg_kr[CLASSES]
    return kr*mask

# def calculate_decay_twigs(input_dst,landmask):
#     print('#########################')
#     print('[[kt ==> DECAY CONSTANT OF LITTER TWIGS]]')
#     print('kt = 0.276')
#     print('Units = yr-1')
#     print('#########################\n')    
#     kt = 0.276 * input_dst * landmask
#     return kt

# def calculate_decay_wood(input_dst,landmask):
#     print('#########################')
#     print('[[kw ==> DECAY CONSTANT OF DEAD WOOD]]')
#     print('kw = 0.083')
#     print('Units = yr-1')
#     print('#########################\n')        
#     kw = 0.083 * input_dst * landmask
#     return kw

# def calculate_litter_production(input_dst,landmask):
#     print('#########################')
#     print('[[lp ==> LITTER PRODUCTION ]]')
#     print('lp = 0.13')
#     print('Units = g/cm2/yr')
#     print('#########################\n')    
#     litter = 0.13 * input_dst * landmask
#     return litter

def get_mineral_deposition(depos,landmask):
    print('#########################')
    print('[[si ==> ANNUAL DEPOSITION OF MINERAL SEDIMENT ON THE SURFACE]]')
    print('si = model deposition * volumetric mass')
    print('Units = g/cm2/yr')
    #from model = cm/yr
    #si = g/cm2/yr
    volumetric_mass = 2650 #kg/m3
    
    #sand_density = 1.52 #g/cm3)
    print('Model ouptut for deposition is cm/year)')
    print('Assuming volumetric mass of %s km/m3' %(volumetric_mass))
    print('Converting to g/cm2/year')
    print('#########################\n')    
    si = depos*volumetric_mass*0.001 
    
    return si*mask
    

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


def get_input_files(AOI,aoi_dir,ref_dir,EPSG,bounds_4326):
	
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
    tmp_dir = aoi_dir / 'TMP'
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)
    ##########################################################################################
    ## AVIRIS-NG Aboveground Biomass maps
    ## Units are Mg/Ha
    print('')
    print('##### ABOVEGROUND BIOMASS\n\n')
    final_file = aoi_dir / 'AVIRIS_AGB.tif'
    
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
    B = np.where( SALINITY == 3, 1,np.where((SALINITY==4) | (SALINITY ==5) , 2, np.where(SALINITY==6,3,np.nan)))
    CLASSES = ((A*3)+B).astype(int)

    ##########################################################################################
    ## mineral deposition results from hydro/morphodynamic models
    ## units cm/year?
    print('')
    print('##### Delft3D Inorganice Mass Accumulation Rate (IMAR)\n\n')
    final_file = aoi_dir / 'Delft3D_IMAR.tif'
   

    print('# Search EarthData for Delft3D Sediment Model Data over the AOI\n')
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
    print('# Crop the Delft3D IMAR file to AOI and resample to %sm \n' %(res))
    os.system('gdalwarp -overwrite -tr %s %s %s/%s_%s.tif %s '\
            ' -te %s %s %s %s -srcnodata -9999 -dstnodata -9999 -co COMPRESS=DEFLATE -q'
            %(res,res,tmp_dir, filename, EPSG,final_file,ulx,lry,lrx,uly))
    print('--> Saved as %s\n' %(final_file))
    print('')


    IMAR = rasterio.open(final_file).read(1)
    ##########################################################################################

    ##########################################################################################
    ## new water mask
    ## 1 = water
    ## 0 = not water
    print('##### Delta-X Watermask')

    watermask_file = ref_dir  / 'Delta-X_watermask.tif'
    os.system('gdalwarp -q -t_srs epsg:%s %s %s' %(EPSG,watermask_file,str(watermask_file)[:-4] + '_%s.tif' %(EPSG)))
    os.system('gdalwarp -overwrite -tr %s %s %s %s '\
            ' -te %s %s %s %s -srcnodata -9999 -dstnodata -9999 -co COMPRESS=DEFLATE -q'
            %(res,res,str(watermask_file)[:-4] + '_%s.tif' %(EPSG),aoi_dir/'watermask.tif',ulx,lry,lrx,uly))
    WATER = rasterio.open(aoi_dir/ 'watermask.tif').read(1)
    watermask = np.where(WATER ==0,1,np.nan)
    ##########################################################################################
    
    mask = np.where((watermask == 1) & (CLASSES >0),1,np.nan)


    return AGB, SALINITY, LANDCOVER, IMAR, BASINS, CLASSES, mask, lats, lons


    
    
     
    
    
def calculate_other_parameters(AVIRIS, SALINITY, LANDCOVER,b):
    ## Aboveground Organics
    # litter
    
    ref_arr = np.ones((AVIRIS.shape))
    c0 = calculate_lignin_leaflitter(ref_arr)
    c1 = calculate_ash_litter(ref_arr)
    f1 = calculate_ratio_leaflitter(ref_arr)
    lp = calculate_litter_production(ref_arr)
    # decay
    ka = calculate_decay_labileOM_surface(ref_arr)
    kt = calculate_decay_twigs(ref_arr)
    kw = calculate_decay_wood(ref_arr)
    ke = calculate_litter_export(ref_arr)
    km = calculate_turnover_wood(ref_arr)    
    
    
    ## Belowground Organics
    r0 = calculate_belowground_biomass(b)
    kb = calculate_decay_labileOM_below(ref_arr)
    kr = calculate_turnover_fineroots(ref_arr)
    
    b0 = calculate_bulkdensity_OM(ref_arr)
    
    
    e = calculate_e_distribution(ref_arr)
    
    f2 = calculate_microbial_resp_decomp(ref_arr)
    f3 = calculate_labile2refractory_decomp(ref_arr)
    
    fc1 = calculate_refractory_fineroot(ref_arr)
    fc2 = calculate_refractory_largeroot(ref_arr)
    
    kc = calculate_decay_refractoryOM(ref_arr)

    return c0, c1, f1, lp, ka, kt, kw, ke, km, r0, kb, kr, b0, e, f2, f3, fc1, fc2, kc
    


    
# working_directory = Path('/Users/Alchrist/Downloads/')
# ref_file = working_directory / 'aoi.tif'
# AOI = rasterio.open(ref_file)
# profile = ref.profile
# transform = profile['transform']
# width = profile['width']
# height = profile['height']
# ulx = transform[2]
# uly = transform[5]
# lrx = ulx + width*transform[0]
# lry = uly + height*transform[4]
# res = profile['transform'][0]
# ref_array = ref.read(1)   


# get_input_files(AOI)














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


