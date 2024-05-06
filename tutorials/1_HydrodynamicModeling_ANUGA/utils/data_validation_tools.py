import pandas as pd
import pytz
from datetime import datetime
from tqdm import notebook
import numpy as np

def process_CRMS_datetimes(date, time, time_zones):
    datetime_series = date + ' ' + time + ' ' + time_zones
    # Define timezone mapping
    timezone_map = {
        'CST': pytz.timezone('America/Chicago'),  # Central Standard Time
        'CDT': pytz.timezone('America/Chicago'),  # Central Daylight Time
    }
    # Convert the series
    for index, value in notebook.tqdm(datetime_series.items(), desc = 'Processing CRMS datetimes', leave=False, total = np.shape(date)[0]):
        # Split the datetime string and the timezone part
        datetime_str, tz_str = value.rsplit(' ', 1)
        # Parse the datetime string
        dt = datetime.strptime(datetime_str, '%m/%d/%Y %H:%M:%S')
        # Localize the datetime with the appropriate timezone
        dt = timezone_map[tz_str].localize(dt)
        # Convert to UTC
        dt_utc = dt.astimezone(pytz.utc)
        # Update the series with UTC datetime
        datetime_series.at[index] = dt_utc
    return datetime_series
