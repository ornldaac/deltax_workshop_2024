# Anuga Hydrodynamic Modeling: Hands On Exercises

### Author information:

**Antoine Soloy, Ph.D.**  
Division 334F, Caltech - Jet Propulsion Laboratory  
4800 Oak Grove Drive, Pasadena, CA, USA 91109-8099  
Contact: antoine.soloy@jpl.nasa.gov  

### Notebook Series Description:

This series of notebooks is an introduction to the implementation of hydrodynamic models to simulate the flow of water in complex vegetated estuarine and deltaic environments. The ANUGA hydrodynamic model is used for demonstration purpose, and applied to a user-defined subregion of interest in the Mississippi River Delta, LA, USA.  

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ornldaac/deltax_workshop_2024/blob/main/tutorials/1_HydrodynamicModeling_ANUGA/[1]ANUGA_DEM_processing.ipynb)[1]ANUGA_DEM_preprocessing.ipynb : Modifications to the Digital Elevation Model to optimize hydrological connectivity.
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ornldaac/deltax_workshop_2024/blob/main/tutorials/1_HydrodynamicModeling_ANUGA/[2]ANUGA_friction_map_generation.ipynb)[2]ANUGA_friction_map_generation.ipynb : Generation of a land & water classification map for friction calibration purpose.  
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ornldaac/deltax_workshop_2024/blob/main/tutorials/1_HydrodynamicModeling_ANUGA/[3]ANUGA_model_run.ipynb)[3]ANUGA_model_run.ipynb : Setting up of the hydrodynamic model and run.  
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ornldaac/deltax_workshop_2024/blob/main/tutorials/1_HydrodynamicModeling_ANUGA/[4]ANUGA_output_analysis.ipynb)[4]ANUGA_output_analysis.ipynb : Output reading and generation of animations.  
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ornldaac/deltax_workshop_2024/blob/main/tutorials/1_HydrodynamicModeling_ANUGA/[5]ANUGA_model_validation.ipynb)[5]ANUGA_model_validation.ipynb : Comparison between in-situ/remote sensing data and model predictions.  

### Necessary Datasets:

Please ensure that all the following datasets were downloaded and placed into the "data" subfolder. Zip files should be unziped in the same folder, keeping the original tree structure intact.

- #### Delta-X Digital Elevation Model (DEM) and Water Mask
    - #### Download link: https://doi.org/10.3334/ORNLDAAC/2181
    - #### Citation: Christensen, A.L., M.W. Denbina, and M. Simard. 2023. Delta-X: Digital Elevation Model, MRD, LA, USA, 2021. ORNL DAAC, Oak Ridge, Tennessee, USA. https://doi.org/10.3334/ORNLDAAC/2181

- #### Delta-X Vegetation Classification Map
    - #### Download link: https://drive.google.com/file/d/1OSyAWekAuXGFz3yLRz61PlLuQgj2FkGa/view?usp=drive_link
    - #### Citation: Jensen, D.J., E. CastaÃƒÂ±eda-Moya, E. Solohin, D.R. Thompson, and M. Simard. 2024. Delta-X AVIRIS-NG L3 Derived Vegetation Types, MRD, Louisiana, USA. ORNL DAAC, Oak Ridge, Tennessee, USA. https://doi.org/10.3334/ORNLDAAC/2352

- #### Delta-X AirSWOT L3 Water Surface Elevations
    - #### Download link: https://doi.org/10.3334/ORNLDAAC/2349
    - #### Citation: Denbina, M.W., M. Simard, and E. Rodriguez. 2024. Delta-X: AirSWOT L3 Water Surface Elevations, MRD, Louisiana, 2021, V2. ORNL DAAC, Oak Ridge, Tennessee, USA. https://doi.org/10.3334/ORNLDAAC/2349

- #### JPL Mississippi River Centerlines
    - #### Download link: https://landscape.jpl.nasa.gov/cgi-bin/data-search.pl
    - #### Citation: Christensen, A.L., Soloy, A., Savelli, R., Moritz, J.M., & Simard, M. (2023b). Centerlines of the Mississippi River (V1.0) [Data file]. Retrieved from https://landscape.jpl.nasa.gov/cgi-bin/data-search.pl

- #### CRMS water level gauges
    - #### Download link: https://drive.google.com/file/d/1WjOSv8XcYZwkN_-oPdsdXH_5Lg2zcp5x/view?usp=drive_link (workshop sample)
    - #### Citation: Coastal Protection and Restoration Authority (CPRA) of Louisiana. 2024. Coastwide Reference Monitoring System-Wetlands Monitoring Data. Retrieved from Coastal Information Management System (CIMS) database. http://cims.coastal.louisiana.gov. Accessed 24 January 2024
