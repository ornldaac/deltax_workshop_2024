import numpy as np
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess
from osgeo import gdal

def GenerateTideGauge(filename,
                      t_start,
                      t_end,
                      t_step,
                      offset = 0,
                      smoothing = False,
                      smoothing_span = .1,
                      hot_start = False,
                      last_frame = 0): # timestep in seconds
    """
    Function to generate a tidal BC from gauge data.
    Specify inputs in settings.py. Returns a function of time.
    """
    # Measured water levels
    tides = pd.read_csv(filename, header=0, names=['datetime', 'WL'])
    tides['datetime'] = pd.to_datetime(tides['datetime'], format="%Y-%m-%d %H:%M:%S", utc=True)
    tides['WL'] = tides['WL'].copy()+offset
    
    duration = (t_end-t_start)
    t_step = t_step

    if hot_start==True:
        t_start = t_start + (t_step*last_frame)
    
    tides['seconds'] = (tides['datetime']-t_start).dt.total_seconds()
    
    if smoothing == True:
        # Smoothen raw discharge
        tides['WL'] = lowess(tides['WL'],  tides['datetime'], is_sorted=True, frac=smoothing_span, it=1)[:,1]


    
    #Set the anuga tide boundary function
    fBC_tides = lambda t: np.interp(t, 
                                    tides['seconds'], 
                                    tides['WL'])
    
    return fBC_tides

def GenerateHydrograph(filename,
                       t_start, 
                       t_end,
                       t_step,
                       offset = 0,
                       smoothing = False,
                       smoothing_span = .25,
                       progressive = True,
                       hot_start = False,
                       last_frame = 0):
    """
    Function to generate a hydrograph from USGS gauge data.
    Specify inputs in settings.py. Returns a function of time.
    When calling, filename either name_WLO, name_ATC, name_FRA
    """
    Q = pd.read_csv(filename, header=0, names=['datetime','Q'])
    Q['datetime'] = pd.to_datetime(Q['datetime'], format="%Y-%m-%d %H:%M:%S%z")

    Q['Q'] = Q['Q']+offset
    
    duration = (t_end-t_start)
    t_step = t_step

    if hot_start==True:
        t_start = t_start + (t_step*last_frame)
    
    Q['seconds'] = (Q['datetime']-t_start).dt.total_seconds()

    if smoothing == True:
        # Smoothen raw discharge
        Q['Q'] = lowess(Q['Q'],  Q['datetime'], is_sorted=True, frac=smoothing_span, it=1)[:,1]

    # Apply progressive discharge on the first 100 timesteps
    if progressive == True:
        p = [Q['Q'][i]*((i+1)/100) for i in range(0,100)]
        Q.loc[0:99,'Q'] = p
        
    #Set the anuga discharge inlet function
    fQ = lambda t: np.interp(t, 
                             Q['seconds'], 
                             Q['Q'])
    
    return fQ

def GenerateTideCosine(amplitude=0.25, period=7.2722e-5, phase=0, offset=0.26):
    """
    Function to generate an artificial tidal BC from a sinusoid. Default is
    cosine with same tidal range as LAWMA, freq of N2.
    Returns a function of time.
    """
    # Default is pseudo tide with same tidal range as LAWMA, freq of N2. Time in sec
    def fBC_tides(t):
        return amplitude*np.cos(t*period + phase) + offset
    return fBC_tides

# ------------------------------------------------------------------------------
# Domain settings functions
# ------------------------------------------------------------------------------
def Raster2Mesh(meshX, meshY, raster_filename):
    """
    Function to grab raster value at each mesh cell centroid.
    Used for friction assignment. Specify raster in settings.py.
    Returns a numpy.ndarray
    """
    src = gdal.Open(raster_filename)
    data_array = np.array(src.GetRasterBand(1).ReadAsArray())
    # Get geographic info
    transform = src.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = -transform[5]

    # Loop and grab all values at mesh coords
    meshVal = np.zeros(len(meshX), dtype=int)
    for ii in list(range(0, len(meshX))):
        col = int((meshX[ii] - xOrigin) / pixelWidth)
        row = int((yOrigin - meshY[ii] ) / pixelHeight)
        try:
            meshVal[ii] = data_array[row][col]
        except IndexError:
            meshVal[ii] = 0
    # Return values
    return meshVal

def AssignFricValue(FricVal, n_array, m_array, h_array, D_array):
    """
    Using Friction ID (output of Raster2Mesh), assign friction parameters.
    All Mannings (n) and Baptist (m, hv, D) values returned.
    Specify coefficients for each class in settings.py
    """
    FricVal = FricVal.astype(int)-1 # Subtract 1 assuming map ID's start at 1
    n = np.zeros_like(FricVal, dtype=float)
    m = np.zeros_like(FricVal, dtype=float)
    hv = np.zeros_like(FricVal, dtype=float)
    D = np.zeros_like(FricVal, dtype=float)
    
    for ii in list(range(len(FricVal))):
        n[ii] = n_array[FricVal[ii]]
        m[ii] = m_array[FricVal[ii]]
        hv[ii] = h_array[FricVal[ii]]
        D[ii] = D_array[FricVal[ii]]
    return n, m, hv, D