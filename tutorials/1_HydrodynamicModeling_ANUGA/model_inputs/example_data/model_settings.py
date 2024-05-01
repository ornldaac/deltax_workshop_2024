# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------
import pandas as pd
import numpy as np
import os
import utm

# Define the path to scripts and data
workshop_dir = os.getcwd()
data_dir = os.path.join(workshop_dir, 'data')
model_inputs_dir = os.path.join(workshop_dir, 'model_inputs')
model_outputs_dir = os.path.join(workshop_dir, 'model_outputs')
model_visuals_dir = os.path.join(workshop_dir, 'visuals')
model_validation_dir = os.path.join(workshop_dir, 'validation')

# Simulation time definition
sim_starttime = pd.to_datetime('2021-03-27 00:00:00', format="%Y-%m-%d %H:%M:%S", utc=True)
sim_endtime = pd.to_datetime('2021-03-28 00:00:00', format="%Y-%m-%d %H:%M:%S", utc=True)
sim_timestep = pd.to_timedelta(5, 'm')
sim_total_duration = (sim_endtime-sim_starttime)
data_download_dates = (str(sim_starttime-pd.to_timedelta(1, 'd'))[0:10].replace(' ', '').replace(':', ''),
                       str(sim_endtime+pd.to_timedelta(1, 'd'))[0:10].replace(' ', '').replace(':', ''))
sim_time = np.arange(sim_starttime, sim_endtime+sim_timestep, sim_timestep)

sim_starttime_str = str(sim_starttime)[0:19].replace('-', '').replace(' ', '').replace(':', '')

# Output file naming
domain_name = 'DX_Workshop'
model_name = f'{sim_starttime_str}_{domain_name}_{sim_total_duration.days}_days'

# Boundary Conditions
discharge_gauge_x, discharge_gauge_y = utm.from_latlon(29.692611, -91.211833)[0:2] # Morgan City
discharge_gauge_ID = ('07381600', 'Morgan City', discharge_gauge_x, discharge_gauge_y)

tide_gauge_x, tide_gauge_y = utm.from_latlon(29.373057, -91.383868)[0:2] # Eugene Island
tide_gauge_ID = ('8764314', 'Eugene Island', tide_gauge_x, tide_gauge_y)

# Gauges filepath
f_discharge =  os.path.join(model_inputs_dir, 'Discharge_at_%s_%s-%s.csv' % (
               discharge_gauge_ID[1],
               data_download_dates[0].replace('-', ''), 
               data_download_dates[1].replace('-', '')))

f_tides =  os.path.join(model_inputs_dir, 'Tides_at_%s_%s-%s.csv' % (
              tide_gauge_ID[1],
              data_download_dates[0].replace('-', ''), 
              data_download_dates[1].replace('-', '')))
