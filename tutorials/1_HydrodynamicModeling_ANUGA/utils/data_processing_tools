import pandas as pd
import geopandas as gpd
import rasterio as rio
import numpy as np
import os
from osgeo import ogr, gdal, osr
from scipy import ndimage
from skimage.measure import label, regionprops_table
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from rasterio.warp import reproject, Resampling, calculate_default_transform
from shapely.ops import linemerge, split
from shapely.geometry import MultiLineString, LineString, MultiPoint, Point, Polygon
from shapely import wkt, geometry, segmentize
from tqdm import notebook, tqdm
import math
import ipyleaflet as ileaf
import ipywidgets as ipyw
from ipyleaflet import Map, DrawControl, Polyline, Polygon as LeafletPolygon
from IPython.display import display, clear_output
import traceback

def setup_line_draw_control():
    """Setup and return a draw control for lines with style options."""
    draw_control = DrawControl()
    draw_control.polyline = {
        "shapeOptions": {"color": "#6b8e23", "weight": 4}
    }
    draw_control.polygon = {}
    draw_control.rectangle = {}
    draw_control.circle = {}
    draw_control.circlemarker = {}
    return draw_control

def save_line_intersection_to_shapefile(filename, line_geometry, polygon_geometry, target_epsg):
    """Calculate the intersection and save the line segment within the polygon to a shapefile."""
    try:
        if line_geometry.is_empty or polygon_geometry.is_empty:
            print("No line or polygon to process. Please draw both.")
            return
        
        # Check intersection
        intersection = line_geometry.intersection(polygon_geometry)
        
        if intersection.is_empty:
            print("The drawn line does not intersect with the polygon.")
            return
        
        # Create a GeoDataFrame for the intersected segment
        gdf = gpd.GeoDataFrame([{'id': 1}], geometry=[intersection], crs="EPSG:4326")
        gdf = gdf.to_crs(epsg=target_epsg)
        gdf.to_file(filename, driver='ESRI Shapefile')
        print(f"Line segment inside the polygon saved to '{filename}' in EPSG:{target_epsg}.")
    except Exception as e:
        print("Failed to process the line and polygon:")
        traceback.print_exc()

def draw_line_and_polygon_interaction(filename="line_within_polygon.shp", epsg=32615, domain_polygon=None):
    center = (29.5, -91.3)
    zoom = 8
    
    map_view = create_map(center, zoom)
    draw_control = setup_line_draw_control()
    map_view.add_control(draw_control)
    
    display_existing_polygon(map_view, domain_polygon)
    
    save_button = ipyw.Button(description="Save Line Inside Polygon")
    output = ipyw.Output()

    polygon_shape = Polygon(domain_polygon) if domain_polygon is not None else None

    def on_button_clicked(b):
        with output:
            clear_output(wait=True)
            if draw_control.last_draw and polygon_shape:
                data = draw_control.last_draw['geometry']['coordinates']
                line = LineString(data)
                save_line_intersection_to_shapefile(filename, line, polygon_shape, epsg)
            else:
                print("Please draw a line and ensure a polygon is provided.")

    save_button.on_click(on_button_clicked)
    
    display(map_view, save_button, output)

def create_map(center, zoom):
    """Create and return a map centered on the given coordinates."""
    return ileaf.Map(center=center, zoom=zoom)

def setup_draw_control():
    """Setup and return a draw control for polygons with style options."""
    draw_control = DrawControl()
    draw_control.polygon = {"shapeOptions": {"color": "#fca103"}}
    return draw_control

def display_existing_polygon(map_view, domain_polygon):
    """Display an existing polygon on the map."""
    if domain_polygon is not None:
        # Ensure the domain_polygon is a list of lists (not a NumPy array or other format)
        domain_polygon = domain_polygon.tolist() if isinstance(domain_polygon, np.ndarray) else domain_polygon
        polygon = LeafletPolygon(
            locations=[{"lat": lat, "lng": lon} for lon, lat in domain_polygon],
            color="blue",
            fill_color="blue"
        )
        map_view.add_layer(polygon)

def save_to_shapefile(filename, geometry, target_epsg):
    """Save polygon geometry to a shapefile with coordinate system transformation."""
    try:
        if geometry.is_empty:
            print("No polygon to save. Please draw a polygon first.")
            return
        
        gdf = gpd.GeoDataFrame([{'id': 1}], geometry=[geometry], crs="EPSG:4326")
        gdf = gdf.to_crs(epsg=target_epsg)
        gdf.to_file(filename, driver='ESRI Shapefile')
        print(f"Polygon coordinates successfully saved to '{filename}' in EPSG:{target_epsg}.")
    except Exception as e:
        print("Failed to save the polygon coordinates:")
        traceback.print_exc()

def draw_polygon(filename="dxws_domain.shp", epsg=32615, domain_polygon=None):
    center = (29.5, -91.3)
    zoom = 8
    
    map_view = create_map(center, zoom)
    draw_control = setup_draw_control()
    map_view.add_control(draw_control)
    
    display_existing_polygon(map_view, domain_polygon)
    
    save_button = ipyw.Button(description="Save to Shapefile")
    output = ipyw.Output()

    def on_button_clicked(b):
        with output:
            clear_output(wait=True)
            if draw_control.last_draw:
                data = draw_control.last_draw['geometry']['coordinates'][0]
                polygon = Polygon(data)
                save_to_shapefile(filename, polygon, epsg)
            else:
                print("Please draw a polygon before saving.")

    save_button.on_click(on_button_clicked)
    
    display(map_view, save_button, output)

def createBuffer(inputfn, outputBufferfn, bufferDist):
    inputds = ogr.Open(inputfn)
    inputlyr = inputds.GetLayer()

    shpdriver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(outputBufferfn):
        shpdriver.DeleteDataSource(outputBufferfn)
    outputBufferds = shpdriver.CreateDataSource(outputBufferfn)
    bufferlyr = outputBufferds.CreateLayer(outputBufferfn, geom_type=ogr.wkbPolygon)
    featureDefn = bufferlyr.GetLayerDefn()

    for feature in inputlyr:
        ingeom = feature.GetGeometryRef()
        geomBuffer = ingeom.Buffer(bufferDist)

        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(geomBuffer)
        bufferlyr.CreateFeature(outFeature)
        outFeature = None

def reproj_match(infile, match, outfile):
    """Reproject a file to match the shape and projection of existing raster. 
    
    Parameters
    ----------
    infile : (string) path to input file to reproject
    match : (string) path to raster with desired shape and projection 
    outfile : (string) path to output file tif
    """
    # open input
    with rio.open(infile) as src:
        src_transform = src.transform
        
        # open input to match
        with rio.open(match) as match:
            dst_crs = match.crs
            
            # calculate the output transform matrix
            dst_transform, dst_width, dst_height = calculate_default_transform(
                src.crs,     # input CRS
                dst_crs,     # output CRS
                match.width,   # input width
                match.height,  # input height 
                *match.bounds,  # unpacks input outer boundaries (left, bottom, right, top)
            )

        # set properties for output
        dst_kwargs = src.meta.copy()
        dst_kwargs.update({"crs": dst_crs,
                           "transform": dst_transform,
                           "width": dst_width,
                           "height": dst_height,
                           "nodata": -9999})
#         print("Coregistered to shape:", dst_height,dst_width,'\n Affine',dst_transform)
        # open output
        with rio.open(outfile, "w", **dst_kwargs) as dst:
            # iterate through bands and write using reproject function
            for i in range(1, src.count + 1):
                reproject(
                    source=rio.band(src, i),
                    destination=rio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=dst_transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.nearest)

def redistribute_vertices(geom, distance):
    if geom.geom_type == 'LineString':
        num_vert = int(round(geom.length / distance))
        if num_vert == 0:
            num_vert = 1
        return LineString(
            [geom.interpolate(float(n) / num_vert, normalized=True)
             for n in range(num_vert + 1)])
    elif geom.geom_type == 'MultiLineString':
        parts = [redistribute_vertices(part, distance)
                 for part in geom]
        return type(geom)([p for p in parts if not p.is_empty])
    else:
        raise ValueError('unhandled geometry %s', (geom.geom_type,))
    return

def resample_gdf(gdf, d=10):
    gdf = gdf.dissolve().reset_index(drop=True)
    gdf = gdf.explode(ignore_index=True).reset_index(drop=True)

    lines = linemerge(gdf['geometry'].to_numpy())
    gdf = gpd.GeoDataFrame(geometry=list(lines.geoms))

    gdf['geometry'] = gdf['geometry'].apply(lambda x: redistribute_vertices(x, d))

    # Segmentize the line strings
    gdf['n_vertices'] = gdf.geometry.apply(lambda x: len(x.coords))
    gdf_to_append = gpd.GeoDataFrame()

    inds = np.where(gdf['n_vertices'] > 2)[0]
    
    N_seg = np.floor(inds.shape[0] / 10000) + 1
    Ni_seg = 1
    
    while len(inds) > 0:
        new_data = []
        for i in inds:
            new_rows = [LineString([gdf['geometry'].iloc[i].coords[j], gdf['geometry'].iloc[i].coords[j+1]])
                        for j in range(len(gdf['geometry'].iloc[i].coords) - 1)]
            new_data.extend(new_rows)
        
        gdf_to_append = pd.concat([gdf_to_append, gpd.GeoDataFrame(geometry=new_data)], ignore_index=True)
        
        gdf = pd.concat([gdf.drop(inds), gdf_to_append], ignore_index=True)
        gdf['n_vertices'] = gdf.geometry.apply(lambda x: len(x.coords))
        gdf_to_append = gpd.GeoDataFrame()

        inds = np.where(gdf['n_vertices'] > 2)[0]
        if len(inds) >= 10000:
            inds = inds[:10000]

    # Remove duplicates
    gdf["reversed"] = gdf["geometry"].apply(lambda x: LineString(x.coords[::-1]))
    gdf = gdf.drop_duplicates(subset=["geometry", "reversed"])
    gdf = gdf.drop(columns=["reversed", "n_vertices"]).reset_index(drop=True)
    
    return gdf

def apply_buffer(BWimage, BWdist, buffer=5):
    mask = BWdist>=buffer
    return mask

def apply_buffer_gen(BWimage, buffer=5):
    BWdist = ndimage.distance_transform_edt(BWimage)
    mask = BWdist>=buffer
    return mask

def BWinvert(BW):
    return (BW-1)*-1

def find_ocean(BWimage, BWdist, resolution_m, buffer_m=1000):
    buffer = np.round(buffer_m/resolution_m)
    large_water_bodies = apply_buffer(BWimage, BWdist, buffer=buffer)
    label_im = label(large_water_bodies)
    properties = ['area']
    labels_table = pd.DataFrame(regionprops_table(label_im, large_water_bodies, 
                  properties=properties))
    ocean_ID = np.where(labels_table['area']==np.max(labels_table['area']))[0][0]+1
    ocean_shrunk = ndimage.binary_fill_holes(1*(label_im==ocean_ID))
    ocean = ndimage.binary_fill_holes(BWinvert(apply_buffer_gen(BWinvert(ocean_shrunk), buffer=buffer)))
    return ocean

def find_channels(BWimage, BWdist, ocean, resolution_m, ths_m, excl_elong = 300):
    ths_px = round(ths_m/resolution_m)
    channels = 1*(BWdist>=ths_px) 
    channels = apply_buffer_gen(BWinvert(channels), ths_px)*1
    channels = BWinvert(channels) - ocean

    label_im = label(channels)
    properties = ['area']
    labels_table = pd.DataFrame(regionprops_table(label_im, channels, 
                                                  properties=properties))
    area_ths = excl_elong * (ths_m ** 2)
    channels_IDs = np.where(labels_table['area']>area_ths)
    channels = ((np.sum(label_im==channels_IDs[0][i]+1 for i in range(channels_IDs[0].shape[0])))>=1)*1
    return channels

def find_channels_range(BWimage, BWdist, ocean, resolution_m, min_ths_m, max_ths_m):
    min_ths_px = round(min_ths_m/resolution_m)
    max_ths_px = round(max_ths_m/resolution_m)
    channels = 1*np.logical_and(BWdist>=min_ths_px, BWdist<=max_ths_m) 
    channels = apply_buffer_gen(BWinvert(channels), min_ths_px)*1
    channels[BWimage==0] = 0
    return channels

def calculate_area_of_equilateral_triangle(side_length):
    """Calculate the area of an equilateral triangle given the side length."""
    return (math.sqrt(3) / 4) * side_length ** 2

def classify_points(domain, roi):
    """
    Classifies points in the domain polygon as either inside or outside the ROI polygon.

    Parameters:
    - domain: Nx2 array-like of x and y coordinates for the domain polygon.
    - roi: Nx2 array-like of x and y coordinates for the region of interest (ROI) polygon.

    Returns:
    - inside_indices: List of indices of domain points that are inside the ROI.
    - outside_indices: List of indices of domain points that are outside the ROI.
    """
    # Create Polygon objects from the points
    domain_polygon = Polygon(domain)
    roi_polygon = Polygon(roi)

    # Initialize lists to hold indices of points inside and outside the ROI
    inside_indices = []
    outside_indices = []

    # Check each point in the domain polygon to see if it is inside the ROI
    for i, point in enumerate(domain_polygon.exterior.coords):
        point_obj = Point(point)
        if roi_polygon.contains(point_obj):
            inside_indices.append(i)
        else:
            outside_indices.append(i)

    return inside_indices, outside_indices
