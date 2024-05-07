import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
import dorado.routines
import rasterio as rio
from tqdm.notebook import tqdm
import dorado.routines as drt
import os

def do_plotting(i, params, walk_data, seed_locs, target_times, filepath):

    # get walk data from most recent time-step
    xia, yia, tia = drt.get_state(walk_data) 
    
    # Initialize figure
    fig = plt.figure(dpi=300, figsize=(5,5))
    ax = fig.add_subplot(111)
    
    # Topography as background
    im_background = ax.imshow(params.topography, vmin=-5, vmax=2, cmap='gist_gray')
    
    # Define the colorbar
    cax = fig.add_axes([ax.get_position().x1+0.01,
                        ax.get_position().y0,0.02,
                        ax.get_position().height])
    cbar = plt.colorbar(im_background, cax=cax)
    cbar.set_label('Bed elevation [m]')
    
    # Plot the seed location
    for ii in range(len(seed_locs)):
        ax.scatter(seed_locs[ii][1], seed_locs[ii][0], c='magenta',
                   edgecolors='black', s=10, linewidths=0.5)

    # Plot the new seed location
    newloc_a = ax.scatter(yia, xia, c='red', s=0.2)
    newloc_a.set_offsets(np.array([yia,xia]).T)
    plt.draw()

    # set the title and save the figure
    ax.set_title('Depth at Hour %s' % int(target_times[i]/3600))
    plt.savefig(filepath, bbox_inches='tight')
    plt.close()    
