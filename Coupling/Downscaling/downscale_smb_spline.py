# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 13:17:26 2023

@author: Kamilla Hauknes Sjursen
"""

#%% Specify python environment
# pyenv
# pyenv('Version','/uio/kant/geo-geofag-u1/henninma/issm/trunk-jpl/projects/jost/py39/bin/python','ExecutionMode','OutOfProcess');

#%% Libraries

# Standard libraries
import os

# External libraries
import numpy as np
import xarray as xr
from scipy import stats
# from scipy.interpolate import interp2d
from scipy.interpolate import RectBivariateSpline
from datetime import datetime
print('All modules loaded')

# Standard libraries
# os=py.importlib.import_module('os');
# 
# # External libraries
# np=py.importlib.import_module('numpy');
# xr=py.importlib.import_module('xarray');
# scipy=py.importlib.import_module('scipy');
# datetime=py.importlib.import_module('datetime');

# Internal libraries

#%% Function downscale_mb_grid

def downscale_mb_grid(da_monthly_mb, ds_dem, ds_dem_hr, ds_glfrac, ds_glfrac_hr):
        
    """
    MÅ OPPDATERES!
    
    Downscales mass balance grid from coarse to fine resolution.

    Uses the statistical downscaling method by Noel et al. (2016). Regression parameters 
    are determined for the equation:

    mb_c = a_c + b_c * h_c
    
    Where mb_c is the mass balance in the coarse grid cell, and h_c is the elevation
    in the coarse grid cell. First b_c is determined by regression to each ice-covered
    cell using the current cell and available neighboring cells (from 3-9 cells in total).
    The a_c is determined by using the b_c together with mb_c in the current grid cell.

    When arrays of a_c and b_c are determined, these are bilinearly interpolated
    to the high-resolution grid to get arrays a_h and b_h. The regression parameters
    a_h and b_h are then used together with the elevation, h_h, in the high-resolution
    DEM to calculate the high-resolution mass balance, mb_h.      
    
    Parameters
    ----------
    da_monthly_mb (year_month,Y,X) : xarray.DataArray (float)
        DataArray of gridded monthly surface mass balance for the 
        glacierized area. Coarse resolution mass balance.
        Coordinates:
        year_month : object
            Time index, month & year. 
        Y : float
            Y-coordinate of cell centers.
        X : float
            X-coordinate of cell centers. 
    ds_dem (Y,X) : xarray.Dataset
        Dataset for the catchment/glacier with coarse resolution containing:
        Coordinates:
        Y : float
            Y-coordinate of cell centers.
        X : float
            X-coordinate of cell centers.
        Data variables:
        elevation (Y,X) : xarray.DataArray (float)
            Elevation in each cell in bounding box.
        glacier_fraction (Y,X) : xarray.DataArray (float)
            Fraction of cell inside glacier boundary.
        Attributes:
        res : float
            Resolution of DEM (cellsize).  
    ds_dem_hr (Y,X) : xarray.Dataset
        Dataset for the catchment/glacier with high resolution containing:
        Coordinates:
        Y : float
            Y-coordinate of cell centers.
        X : float
            X-coordinate of cell centers.
        Data variables:
        elevation (Y,X) : xarray.DataArray (float)
            Elevation in each cell in bounding box.
        glacier_fraction (Y,X) : xarray.DataArray (float)
            Fraction of cell inside glacier boundary.
        Attributes:
        res : float
            Resolution of DEM (cellsize).  

    Returns
    -------
    da_monthly_mb_hr (year_month,Y,X) : xarray.DataArray (float)
        DataArray of gridded monthly surface mass balance for the 
        glacierized area. High resolution mass balance.
        Coordinates:
        year_month : object
            Time index, month & year. 
        Y : float
            Y-coordinate of cell centers.
        X : float
            X-coordinate of cell centers. 
    """

    def get_mb_regression(mb, h, gl_mask_1, r_win, c_win):

        """
        Find regression coefficients a and b in the
        equation y = a + bx. First find b by regression to
        current cell and and adjoining cells (3-9 cells in total).
        Then find a by using y and h in current cell and
        a = y - bx.
    
        Parameters
        ----------
        mb : np.array
            Array with values on which to perform regression (y).
        h : np.array
            Array with x variable.
        gl_mask_1: np.array
            Array with glacier mask to determine cells part of
            glacierized area.
    
        Returns
        -------
        a : np.array
            Array of regression coefficient a.
        b : np.array
            Array of regression coefficient b.
        """
        rows, cols = gl_mask_1.shape
        indices = np.arange(0, rows * cols).reshape(rows, cols)
        gl_mask_nonan = gl_mask_1.copy()
        gl_mask_nonan[np.isnan(gl_mask_1)] = 0
        indices[gl_mask_nonan < 1] = 0
        indices_flat = indices[indices != 0]
    
        mb_flat = mb.flatten()
        h_flat = h.flatten()
        b = np.empty((mb_flat.shape))
        b.fill(np.nan)
        a = np.empty((mb_flat.shape))
        a.fill(np.nan)
    
        for i in range(0, len(indices_flat)):
            idx = indices_flat[i]
            indices_sub = np.array([idx-1, idx, idx+1, (idx+cols-1),(idx+cols),(idx+cols+1), (idx-cols-1),(idx-cols),(idx-cols+1)])
            indices_sub = indices_sub[np.isin(indices_sub, indices_flat)]
            mb_sub = mb_flat[indices_sub]
            h_sub = h_flat[indices_sub]
    
            res = stats.linregress(h_sub, mb_sub)
            b[idx] = res.slope
            a[idx] = mb_flat[idx] - b[idx] * h_flat[idx]
    
        b = b.reshape(rows,cols)
        a = a.reshape(rows,cols)
    
        return a, b

    def pad_regression(arr, r_win, c_win):
        rows, cols = arr.shape
            
        for i in range(2):
            for y in range(0, rows):
                
                ymin = max(0, y - r_win)
                ymax = min(rows, y + r_win + 1)
                
                for x in range(0, cols):
                    
                    xmin = max(0, x - c_win)
                    xmax = min(cols, x + c_win + 1)
                    check_mask = arr[y, x]
                    
                    if np.isnan(check_mask):
                        
                        arr_sub = arr[ymin:ymax, xmin:xmax]
                        
                        if np.sum(~np.isnan(arr_sub)) >= 3:
                            
                            arr[y,x] = np.mean(arr_sub[~np.isnan(arr_sub)])
        return arr

    def downscale_coeff(arr, Y_src, X_src, Y_dest, X_dest, scaling):
    
        x_src = np.arange(0, len(X_src), 1)
        y_src = np.arange(0, len(Y_src), 1)
        x_dest = np.arange(0, len(X_src), scaling)
        y_dest = np.arange(0, len(Y_src), scaling)
    
        f = RectBivariateSpline(y_src, x_src, arr)
        arr_interp = f(y_dest, x_dest)
        return(arr_interp)
    
    print('starting_downscaling')
    
    # Get coarse elevation.
    elev = ds_dem.surface_senorge.values    
    #elev[elev==0] = np.nan
    
    # Get coarse glacier mask.
    mask_glacierized_area = np.array(ds_glfrac.dummy_glacier_fraction)
    
    # Get coarse mask of ones for cells part of glacier, zero otherwise.
    mask_glacier_1_gridded = mask_glacierized_area.copy()
    # RuntimeWarning: invalid value encountered in greater:
    mask_glacier_1_gridded[mask_glacier_1_gridded>0]=1
    
    # Get high resolution elevation. 
    elev_h = np.array(ds_dem_hr.surface_modelled.values)
    gl_frac_h = np.array(ds_glfrac_hr.dummy_glacier_fraction.values)
    gl_frac_h[gl_frac_h == 0] = np.nan
    
    # Get source and destination coordinates.
    Y_c = ds_dem.Y.values
    X_c = ds_dem.X.values
    Y_h = ds_dem_hr.Y.values
    X_h = ds_dem_hr.X.values
    
    # Get source and destination cellsize.
    cellsize_c = ds_dem.res
    cellsize_h = ds_dem_hr.res
    coarseness = cellsize_h/cellsize_c
    
    # Set window shape for padding and regression.
    shape_win = (3, 3)
    row_win = np.floor(shape_win[0] / 2).astype(int)
    col_win = np.floor(shape_win[1] / 2).astype(int)
    
    # Create empty array to fill high resolution monthly mass balance.
    monthly_mb_h = np.empty((len(da_monthly_mb.year_month), len(Y_h), len(X_h)))
    monthly_mb_h.fill(np.nan)
    
    for i in range(0, len(da_monthly_mb.year_month)):
        
        # Get coarse resolution gridded monthly mass balance.
        mb_c = da_monthly_mb[i,:,:].values
        
        # We have a coarse resolution elevation grid (Y1,X1). 
        # For array mb_c together with elevation (Y1,X1) we compute the
        # regression slope b_mb_c and intercept a_mb_c in each cell that is glacier covered
        # Arrays contain regression parameters in cells that are glacier-covered,
        # otherwise nan. 
        a_mb_c, b_mb_c = get_mb_regression(mb_c, elev, mask_glacier_1_gridded, row_win, col_win)
        
        # Pad regression parameters 
        # Padded twice.
        a_mb_c = pad_regression(a_mb_c, row_win, col_win)
        b_mb_c = pad_regression(b_mb_c, row_win, col_win)
    
        # Replace Nan with zeros.
        a_mb_c[np.isnan(a_mb_c)] = 0
        b_mb_c[np.isnan(b_mb_c)] = 0
        
        # Downscale coefficients to high-resolution grid.    
        a_mb_h = downscale_coeff(a_mb_c, Y_c, X_c, Y_h, X_h, coarseness)
        b_mb_h = downscale_coeff(b_mb_c, Y_c, X_c, Y_h, X_h, coarseness)
    
        # Calculate high resolution mass balance in each cell.
        mb_h = a_mb_h + b_mb_h * elev_h
    
        # Multiply with glacier mask to mask cells not part of domain with nan.
        mb_h = mb_h * gl_frac_h
        
        # Store in array.
        monthly_mb_h[i,:,:] = mb_h
    
    # Convert to DataArray.
    da_monthly_mb_h = xr.DataArray(monthly_mb_h,
                                 coords= {'year_month': da_monthly_mb.year_month.values,
                                          'Y': Y_h,
                                          'X': X_h},
                                 dims=["year_month", "Y", "X"],
                                 attrs={'Name': 'Gridded monthly mass balance',
                                        'res': cellsize_h,
                                        'dim': 'mm w.e.'})
    
    # Example plot 
    # da_monthly_mb_h[0,:,:].plot(cmap=plt.cm.RdBu, #vmin = -1000, vmax = 1000,
    #                           cbar_kwargs={'label':'Monthly mass balance (mm w.e.)'})
    # plt.xlabel('$UTM_{east}$')
    # plt.ylabel('$UTM_{north}$')
    # plt.title('Jan 1960')
    
    return da_monthly_mb_h

# End of function downscale_mb_grid()

#%% Read files and run downscaling function.

# def run_smb_downscaling(start_yr, end_yr):
def run_smb_downscaling(start_yr, end_yr, path_smb, model_string):
        
    # Filepath and filenames of DEMs and SMB files.
    filepath = '/uio/hypatia/geofag-personlig/geohyd-staff/henninma/issm/trunk-jpl/projects/jost/Coupling/'
    filepath_jost = '/uio/hypatia/geofag-personlig/geohyd-staff/henninma/issm/trunk-jpl/projects/jost/'
    filepath_output_ds = '/uio/hypatia/geofag-personlig/geohyd-staff/henninma/issm/trunk-jpl/projects/jost/Coupling/Downscaling/Output/'
#     filepath = '/uio/kant/geo-geofag-u1/henninma/issm/trunk-jpl/projects/jost/Coupling/'
#     filepath_output_ds = '/uio/kant/geo-geofag-u1/henninma/issm/trunk-jpl/projects/jost/Coupling/OutputDownscaled/Test/'

    # Make new folder in directory with 
#     newdir = os.path.join(filepath_output_ds, datetime.now().strftime('%Y-%m-%d_%H_%M_%S'))
#     os.makedirs(newdir)
#     newdir = filepath_output_ds + 'Test'

    filename_coarse_dem = 'Downscaling/base_files/dem_JOB_ref_SMB_1km_UTM32.nc'
    filename_base_fine_dem = 'Downscaling/base_files/dem1966_JOB_ref_SMB_100m_UTM32.nc'
#     filename_issm_output_dem = 'issmoutput.nc'
    filename_issm_output_dem = ('Output/' + model_string + 'issmoutput.nc')
    
    filename_coarse_glfrac = 'Downscaling/base_files/gl_frac_dummy_JOB_ref_SMB_1km_UTM32.nc'
    filename_fine_glfrac = 'Downscaling/base_files/gl_frac_dummy_JOB_ref_SMB_100m_UTM32.nc'

#     filename_monthly_mb_coarse = 'Downscaling/base_files/monthly_refmb_JOB_1km_1960_2020_UTM32.nc'
    filename_monthly_mb_coarse = (filepath_jost + path_smb)
#     filename_monthly_mb_coarse = 'Downscaling/base_files/monthly_refmb_JOB_1km_1960_2020_UTM32_Jan24_Pcorr_glacieronly.nc'

    # Read DEMs and SMB files as xarray Datasets.
    # Coarse resolution DEM and dummy glacier fraction.
    with xr.open_dataset(filepath + filename_coarse_dem) as ds_coarse_dem_out:
        ds_coarse_dem = ds_coarse_dem_out
 
    ds_coarse_dem.attrs = {'res': 1000}
        
    with xr.open_dataset(filepath + filename_coarse_glfrac) as ds_coarse_glfrac_out:
        ds_coarse_glfrac = ds_coarse_glfrac_out
    
    # High resolution base DEM and dummy glacier fraction
    with xr.open_dataset(filepath + filename_base_fine_dem) as ds_base_fine_dem_out:
        ds_base_fine_dem = ds_base_fine_dem_out
        
    with xr.open_dataset(filepath + filename_fine_glfrac) as ds_fine_glfrac_out:
        ds_fine_glfrac = ds_fine_glfrac_out
    
    # Modelled surface (output from ISSM).
    with xr.open_dataset(filepath + filename_issm_output_dem) as ds_issm_dem_out:
        ds_issm_dem_t = ds_issm_dem_out
        
    ds_issm_dem = ds_issm_dem_t.transpose()    
    ds_issm_dem.attrs = {'res': 100}
    
    # Overlay new DEM of glacier area with values in 1966 DEM
    # Where values in ds_issm_dem are nan, values from ds_base_fine_dem 
    # are inserted.
    ds_fine_dem = ds_issm_dem.combine_first(ds_base_fine_dem)
        
    # Coarse resolution mass balance data.
#     with xr.open_dataset(filepath + filename_monthly_mb_coarse) as ds_mb_coarse_out:
    with xr.open_dataset(filename_monthly_mb_coarse) as ds_mb_coarse_out:
        ds_mb_coarse = ds_mb_coarse_out
        
    # Define year and month from start and end years.
    start_yrmonth = str(start_yr) + '-01'
    end_yrmonth = str(end_yr) + '-12'

    # Get monthly smb from file .
    da_mb_coarse = ds_mb_coarse.mb_monthly.sel(year_month = slice(start_yrmonth, end_yrmonth))

    # Downscale smb for start year to end year from 1km to 100km grid 
    da_monthly_mb_fine = downscale_mb_grid(da_mb_coarse, ds_coarse_dem, ds_fine_dem, ds_coarse_glfrac, ds_fine_glfrac)

    # Convert to Dataset
    ds_monthly_mb_fine = da_monthly_mb_fine.to_dataset(promote_attrs=True, name='mb_monthly')



    # Save modelled surface and modelled mb
    filename_monthly_mb_fine = model_string + 'smb_100m_' + str(start_yr) + '_' + str(end_yr) + '.nc'
    filename_surface_modelled = model_string + 'surface_modelled_100m_' + str(start_yr-1) + '.nc'
    filename_smb_output = model_string + 'smboutput.nc'

    # Save one version of surface elevation, and two versions
    # of mass balance (one for storage and one for ISSM input).
#     ds_fine_dem.to_netcdf(newdir + '/' + filename_surface_modelled)
#     ds_monthly_mb_fine.to_netcdf(newdir + '/' + filename_monthly_mb_fine)
    
    ds_fine_dem.to_netcdf(filepath_output_ds + filename_surface_modelled)
    ds_monthly_mb_fine.to_netcdf(filepath_output_ds + filename_monthly_mb_fine)
    ds_monthly_mb_fine.to_netcdf(filepath + '/Output/' + filename_smb_output)
#     ds_monthly_mb_fine.to_netcdf(filepath + filename_smb_output)
# 
#     with open(filepath_output_ds + filename_surface_modelled, 'wb') as f:  ds_fine_dem.to_netcdf(f)

    finished_downscaling = True

    return finished_downscaling

#%% Run downscaling

# Run manually:
#done = run_smb_downscaling(1966,1975)

# Run in matlab by specifying:
# run = pyrunfile("downscale_smb.py","done", start_year=XXXX, end_year=YYYY)

# run = run_smb_downscaling(start_year, end_year)
run = run_smb_downscaling(start_year, end_year, path_smb, model_string)

#%%
# from matplotlib import pyplot as plt

# # Example plot 
# smb_downscaled[0,:,:].plot(cmap=plt.cm.RdBu, #vmin = -1000, vmax = 1000,                           
#                            cbar_kwargs={'label':'Monthly mass balance (mm w.e.)'})
# plt.xlabel('$UTM_{east}$')
# plt.ylabel('$UTM_{north}$')
# plt.title('Jan 1960')


