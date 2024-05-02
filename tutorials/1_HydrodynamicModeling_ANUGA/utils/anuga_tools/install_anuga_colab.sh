#!/bin/bash

echo "(1) Install pip packages"
pip -q install dill gitpython meshpy netcdf4 Pmw pymetis numpy==1.23 scipy mpi4py scikit-image pandas gmsh-interop triangle rasterio geopandas ipyleaflet ipywidgets cmocean statsmodels dataretrieval noaa_coops utm Cython > /dev/null 2>&1 

echo "(2) Install gdal"
apt-get -q -y install python-gdal gdal-bin  > /dev/null 2>&1 

echo "(3) Install netcdf4"
apt-get -q -y install python-netcdf4  > /dev/null 2>&1 

echo "(4) Install ffmpeg"
apt-get -q -y install python-ffmpeg  > /dev/null 2>&1 

echo "(5) Download anuga_core github repository"
git clone --quiet https://github.com/GeoscienceAustralia/anuga_core.git  > /dev/null 2>&1 

echo "(6) Install anuga"
cd anuga_core
python setup.py --quiet build  > /dev/null 2>&1 
python setup.py --quiet install  > /dev/null 2>&1 
cd ../
rm -rf anuga_core

echo "(7) Download dorado github repository"
git clone --quiet https://github.com/passaH2O/dorado.git  > /dev/null 2>&1 

echo "(8) Install dorado"
cd dorado
python setup.py --quiet build  > /dev/null 2>&1 
python setup.py --quiet install  > /dev/null 2>&1 
cd ../
rm -rf dorado

echo "(9) Ready to go"
