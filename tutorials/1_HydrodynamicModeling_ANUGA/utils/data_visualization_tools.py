import pandas as pd
import numpy as np
from tqdm import notebook, tqdm
from shapely.geometry import Point, Polygon
import matplotlib.tri as tri
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from IPython import display
import rasterio as rio

def euclidian_distance(x1, y1, x2, y2):
    dist = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    return dist

def get_output_values(ds, field, times, locations):
    return ds[field][times, locations].data

def compute_velocity(ds, times, locations):
    stage = get_output_values(ds, 'stage', times, locations)
    elevation = ds['elevation'][locations].data
    height = stage-elevation
    xmom = get_output_values(ds, 'xmomentum', times, locations)
    ymom = get_output_values(ds, 'ymomentum', times, locations)
    
    # Get vel
    minimum_allowed_height=1.0e-03
    h_inv = 1/height
    hWet = (height > minimum_allowed_height)
    xvel = xmom*h_inv*hWet
    yvel = ymom*h_inv*hWet
    vel = (xmom**2 + ymom**2)**0.5*h_inv*hWet
    return xvel, yvel, vel

def convert_seconds(sec):
    minutes, seconds = divmod(int(sec), 60)
    hours, minutes = divmod(int(minutes), 60)
    days, hours = divmod(int(hours), 24)
    return f"{int(days)} days, {int(hours):02}:{int(minutes):02}:{int(seconds):02}"

def animate_stage(saved_output, domain_polygon, gauge_location_ind = None, gauge_time = None, gauge_stage = None, bg_img_filepath = None, fps = 5, initial_frame = 0, final_frame = 1, framestep = 1, dpi = 150, save = False, filename = 'stage.mp4'):
    
    frames = np.arange(initial_frame, final_frame, framestep)
    
    pbar_render = notebook.tqdm(total=len(frames), file=sys.stdout, desc = 'Rendering...')
    
    if save:
        pbar_save = notebook.tqdm(total=len(frames), file=sys.stdout, desc='Saving...')
    
    total_ylim = np.array([np.nanmin(domain_polygon[:,1]), np.nanmax(domain_polygon[:,1])])
    total_xlim = np.array([np.nanmin(domain_polygon[:,0]), np.nanmax(domain_polygon[:,0])])

    x = saved_output['x'][:].data+saved_output.xllcorner
    y = saved_output['y'][:].data+saved_output.yllcorner
    triangles = saved_output['volumes'][:,:].data
    
    Npt = x.shape[0]
    if gauge_time is None:
        gauge_time = np.asarray(saved_output['time'][:].data, dtype=int)
        fig, axes = plt.subplots(2, 1, figsize=(10, 12), dpi=dpi, gridspec_kw={'height_ratios': [100, 0]})
    else:
        fig, axes = plt.subplots(2, 1, figsize=(10, 10), dpi=dpi, gridspec_kw={'height_ratios': [100, 30]})
    
    if bg_img_filepath is not None:
        bg_img_src = rio.open(bg_img_filepath)
        bg_img = np.stack([bg_img_src.read(i) for i in range(1, bg_img_src.count+1)], axis=2)
        bg_img_extent = [bg_img_src.bounds.left, bg_img_src.bounds.right,
                         bg_img_src.bounds.bottom, bg_img_src.bounds.top]
        axes[0].imshow(bg_img, extent = bg_img_extent)
    
    st = get_output_values(saved_output, 'stage', initial_frame, np.arange(0, Npt, 1))
    depth = st - saved_output['elevation'][:].data
    st[depth<=.05] = np.nan

    if gauge_stage is None:
        subfig1 = axes[0].tripcolor(x, y, st, triangles=triangles, cmap='turbo', vmin=-.5, vmax=.5)
        axes[1].axis('off')
    else:
        subfig1 = axes[0].tripcolor(x, y, st, triangles=triangles,
                                   cmap='turbo', vmin=np.nanquantile(gauge_stage, .05), vmax=np.nanquantile(gauge_stage, .95))

    if gauge_location_ind is not None:
        axes[0].scatter(gauge_location_ind[0], gauge_location_ind[1], marker='x', color='w', s=50)
    axes[0].set_aspect('equal', adjustable='box')
    if type(gauge_time[0]) is np.int64 or type(gauge_time[0]) is np.float64:
        t = axes[0].set_title(f'Water level at {convert_seconds(gauge_time[initial_frame])}')
    else:
        t = axes[0].set_title(f'Water level on {str(gauge_time[initial_frame]).replace("T", " ")[:-7]}')
    divider = make_axes_locatable(axes[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(subfig1, cax=cax)
    cbar.set_label('Water level [m]', rotation=90, fontsize=18)
    axes[0].set_xlabel('X', fontsize=18)
    axes[0].set_ylabel('Y', fontsize=18)
    axes[0].set_facecolor('slategrey')
    axes[0].grid(True)
    axes[0].set_xlim(total_xlim)
    axes[0].set_ylim(total_ylim)

    if gauge_stage is not None:
        time_diffs = np.diff(gauge_time)
        if type(gauge_time[0]) is np.int64 or type(gauge_time[0]) is np.float64:
            average_time_step = np.mean(time_diffs) 
            plot_x_increment = average_time_step * framestep
        else:
            average_time_step = np.nanmean(time_diffs.astype('timedelta64[s]')).astype(float)
            # plot_x_increment = average_time_step * framestep
            plot_x_increment = average_time_step*0.003472222208074527/(5*60) #value determined empirically

        subfig2 = axes[1].plot(gauge_time, gauge_stage, 'k-')
        vl = axes[1].vlines(gauge_time[initial_frame], np.min(axes[1].get_ylim()), np.max(axes[1].get_ylim()), colors='r')
        axes[1].set_xlabel('Time')
        axes[1].set_ylabel('Water Level [m]')
        axes[1].tick_params(axis='x', labelrotation = 30)
        axes[1].grid(True)
        v0 = vl.get_segments()[0][0,0]
    
    def updateALL(frame_num):
        st = get_output_values(saved_output, 'stage', frame_num, np.arange(0, Npt, 1))
        depth = st - saved_output['elevation'][:].data
        st[depth<=.05] = np.nan
        subfig1.set_array(np.mean(st[triangles], axis=1))
        if gauge_stage is not None:
            vl.set_segments([np.array([[v0+frame_num*plot_x_increment, np.min(axes[1].get_ylim())], 
                                      [v0+frame_num*plot_x_increment, np.max(axes[1].get_ylim())]])])
        if type(gauge_time[0]) is np.int64 or type(gauge_time[0]) is np.float64:
            t = axes[0].set_title(f'Water level at {convert_seconds(gauge_time[frame_num])} ')
        else:
            t = axes[0].set_title(f'Water level on {str(gauge_time[frame_num]).replace("T", " ")[:-7]}')
        pbar_render.update(1)
        return axes
    
    anim = FuncAnimation(fig, updateALL, frames=frames, interval = 1000 / fps, blit = False, cache_frame_data = False, repeat = True) #
    
    video = anim.to_html5_video()
    html = display.HTML(video)
    display.display(html)
    pbar_render.close()
    plt.close()

    if save:
        FFwriter = animation.FFMpegWriter(fps=fps)
        anim.save(filename, writer = FFwriter, progress_callback = lambda i, n: pbar_save.update(1))
        pbar_save.close()



def animate_velocity(saved_output, domain_polygon, querry_xy_inds, colormap = 'turbo', vcolormin=0, vcolormax=0.5, arrowscale=50, arrowwidth=0.002, arrowheadwidth=5, arrowheadlength=7, arrowheadaxislength=5, gauge_location_ind = None, gauge_time = None, gauge_stage = None, bg_img_filepath = None, fps = 5, initial_frame = 0, final_frame = 1, framestep = 1, dpi = 150, save = False, filename = 'velocity.mp4'):
    
    frames = np.arange(initial_frame, final_frame, framestep)
    
    pbar_render = notebook.tqdm(total=len(frames), file=sys.stdout, desc = 'Rendering...')
    
    if save:
        pbar_save = notebook.tqdm(total=len(frames), file=sys.stdout, desc='Saving...')
    
    total_ylim = np.array([np.nanmin(domain_polygon[:,1]), np.nanmax(domain_polygon[:,1])])
    total_xlim = np.array([np.nanmin(domain_polygon[:,0]), np.nanmax(domain_polygon[:,0])])

    x = saved_output['x'][:].data+saved_output.xllcorner
    y = saved_output['y'][:].data+saved_output.yllcorner
    triangles = saved_output['volumes'][:,:].data
    triang = tri.Triangulation(x, y, triangles)

    Npt = x.shape[0]
    if gauge_time is None:
        gauge_time = np.asarray(saved_output['time'][:].data, dtype=int)
        fig, axes = plt.subplots(2, 1, figsize=(10, 12), dpi=dpi, gridspec_kw={'height_ratios': [100, 0]})
    else:
        fig, axes = plt.subplots(2, 1, figsize=(10, 10), dpi=dpi, gridspec_kw={'height_ratios': [100, 30]})
    
    if bg_img_filepath is not None:
        bg_img_src = rio.open(bg_img_filepath)
        bg_img = np.stack([bg_img_src.read(i) for i in range(1, bg_img_src.count+1)], axis=2)
        bg_img_extent = [bg_img_src.bounds.left, bg_img_src.bounds.right,
                         bg_img_src.bounds.bottom, bg_img_src.bounds.top]
        axes[0].imshow(bg_img, extent = bg_img_extent)
    
    u, v, o = compute_velocity(saved_output, initial_frame, np.arange(0, Npt, 1))
    st = get_output_values(saved_output, 'stage', initial_frame, np.arange(0, Npt, 1))
    depth = st - saved_output['elevation'][:].data
    depth_c = np.nanmean(depth[triangles], axis=1)
    u[depth<=.05] = np.nan
    v[depth<=.05] = np.nan
    o[depth<=.05] = np.nan

    subfig1 = axes[0].tripcolor(triang, np.mean(o[triangles], axis=1), antialiased=True, 
             cmap=colormap, vmin=vcolormin, vmax=vcolormax, alpha=.5, edgecolors='face', linewidths=0)
   
    if gauge_stage is None:
        axes[1].axis('off')

    if gauge_location_ind is not None:
        axes[0].scatter(gauge_location_ind[0], gauge_location_ind[1], marker='x', color='w', s=50)
    axes[0].set_aspect('equal', adjustable='box')
    if type(gauge_time[0]) is np.int64 or type(gauge_time[0]) is np.float64:
        t = axes[0].set_title(f'Water Velocity at {convert_seconds(gauge_time[initial_frame])}')
    else:
        t = axes[0].set_title(f'Water Velocity on {str(gauge_time[initial_frame]).replace("T", " ")[:-7]}')
    divider = make_axes_locatable(axes[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(subfig1, cax=cax)
    cbar.set_label('Water Level [m]', rotation=90, fontsize=18)
    axes[0].set_xlabel('X', fontsize=18)
    axes[0].set_ylabel('Y', fontsize=18)
    axes[0].set_facecolor('slategrey')
    axes[0].grid(True)
    axes[0].set_xlim(total_xlim)
    axes[0].set_ylim(total_ylim)

    # Create the quivers:
    norm = matplotlib.colors.Normalize(vmin=0, vmax=.3)
    cm = matplotlib.cm.CMRmap
    norm_u = u/o
    norm_v = v/o
    noval_inds = np.argwhere((o < 0.01) & (depth < 0.1))
    norm_u[noval_inds] = np.nan
    norm_v[noval_inds] = np.nan

    Quiv = axes[0].quiver(np.nanmean(x[querry_xy_inds], axis=1), 
                          np.nanmean(y[querry_xy_inds], axis=1), 
                          np.nanmean(norm_u[querry_xy_inds], axis=1), 
                          np.nanmean(norm_v[querry_xy_inds], axis=1), 
                          scale=arrowscale, width=arrowwidth, headwidth=arrowheadwidth, headlength=arrowheadlength, headaxislength=arrowheadaxislength,
                          pivot='tail', alpha=.7)

    if gauge_stage is not None:
        time_diffs = np.diff(gauge_time)
        if type(gauge_time[0]) is np.int64 or type(gauge_time[0]) is np.float64:
            average_time_step = np.mean(time_diffs) 
            plot_x_increment = average_time_step * framestep
        else:
            average_time_step = np.nanmean(time_diffs.astype('timedelta64[s]')).astype(float)
            # plot_x_increment = average_time_step * framestep
            plot_x_increment = average_time_step*0.003472222208074527/(5*60) #value determined empirically

        subfig2 = axes[1].plot(gauge_time, gauge_stage, 'k-')
        vl = axes[1].vlines(gauge_time[initial_frame], np.min(axes[1].get_ylim()), np.max(axes[1].get_ylim()), colors='r')
        axes[1].set_xlabel('Time')
        axes[1].set_ylabel('Water Velocity [m/s]')
        axes[1].tick_params(axis='x', labelrotation = 30)
        axes[1].grid(True)
        v0 = vl.get_segments()[0][0,0]
    
    def updateALL(frame_num):
        u, v, o = compute_velocity(saved_output, frame_num, np.arange(0, Npt, 1))
        st = get_output_values(saved_output, 'stage', frame_num, np.arange(0, Npt, 1))
        depth = st - saved_output['elevation'][:].data
        u[depth<=.05] = np.nan
        v[depth<=.05] = np.nan
        o[depth<=.05] = np.nan
        
        subfig1.set_array(np.mean(o[triangles], axis=1))
        
        norm_u = u/o
        norm_v = v/o
        noval_inds = np.argwhere((o < 0.01) & (depth < 0.1))
        norm_u[noval_inds] = np.nan
        norm_v[noval_inds] = np.nan
    
        Quiv.set_UVC(np.nanmean(norm_u[querry_xy_inds], axis=1), np.nanmean(norm_v[querry_xy_inds], axis=1))
        if gauge_stage is not None:
            vl.set_segments([np.array([[v0+frame_num*plot_x_increment, np.min(axes[1].get_ylim())], 
                                      [v0+frame_num*plot_x_increment, np.max(axes[1].get_ylim())]])])
        if type(gauge_time[0]) is np.int64 or type(gauge_time[0]) is np.float64:
            t = axes[0].set_title(f'Water Velocity at {convert_seconds(gauge_time[frame_num])} ')
        else:
            t = axes[0].set_title(f'Water Velocity on {str(gauge_time[frame_num]).replace("T", " ")[:-7]}')
        pbar_render.update(1)
        return axes
    
    anim = FuncAnimation(fig, updateALL, frames=frames, interval = 1000 / fps, blit = False, cache_frame_data = False, repeat = True) #
    
    video = anim.to_html5_video()
    html = display.HTML(video)
    display.display(html)
    pbar_render.close()
    plt.close()

    if save:
        FFwriter = animation.FFMpegWriter(fps=fps)
        anim.save(filename, writer = FFwriter, progress_callback = lambda i, n: pbar_save.update(1))
        pbar_save.close()