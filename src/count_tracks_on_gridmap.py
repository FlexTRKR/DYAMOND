"""
Counts unique tracks on a grid map and save output to a netCDF file.
"""
import numpy as np
import pandas as pd
import xarray as xr
from itertools import chain
import sys, os, time

def get_unique_placements(track_num, lat, lon, ntimes):
    """
    Get unique pairs of lat/lon for a track
    """
    # Put all lat/lon pairs over each time for a track into an array
    this_row = np.array([[lat[track_num, tt], lon[track_num, tt]] for tt in range(0, ntimes)])
    # Return the unique pairs (axis=0)
    return np.array(np.unique(this_row, axis=0))

def count_unique_tracks(lat, lon, xbins, ybins):
    
    ntracks, ntimes = lat.shape

    # A function to loop over tracks
    get_unique = lambda D: get_unique_placements(D, lat, lon, ntimes)

    # Loop over each track and get the unique pairs of lat/lon
    all_uniques = list(map(get_unique, np.arange(0, ntracks)))

    # Flatten the list of lat/lon pairs (using chain), and convert into an array
    unique_latlon = np.array(list(chain(*all_uniques)))

    # Count number on map using histogram2d
    ranges = [[min(ybins), max(ybins)], [min(xbins), max(xbins)]]
    hist2d, yedges, xedges = np.histogram2d(unique_latlon[:,0], unique_latlon[:,1], bins=[ybins, xbins], range=ranges)

    return hist2d


if __name__ == "__main__":

    run_name = sys.argv[1]
    # run_name = 'scream'

    # Sepcify time period
    sdate = '2020-02-01T00'
    edate = '2020-03-01T00'
    # Specify grid
    ranges = [[-60,60], [0,360]]
    xbins = np.arange(0, 360.1, 1)
    ybins = np.arange(-60., 60.1, 1)

    test = ''
    rootdir = '/global/cfs/cdirs/m1867/zfeng/dyamond-winter/'
    indir = f'{rootdir}{run_name}/stats{test}/'
    modfile = f'{indir}trackstats_20200120.0000_20200301.0000.nc'
    outdir = f'{rootdir}{run_name}/stats{test}/monthly/'

    # Convert dates to Pandas datetime
    start_date = pd.to_datetime(sdate)
    end_date = pd.to_datetime(edate)
    # Make output datetime string (yyyymmdd.hhmm)
    sdate_str = start_date.strftime('%Y%m%d')
    edate_str = end_date.strftime('%Y%m%d')

    # Make output filename
    output_filename = f'{outdir}track_counts_gridmap_{sdate_str}_{edate_str}.nc'

    # Read input data
    ds = xr.open_dataset(modfile)
    # Get track initial time valuees
    base_time = ds.base_time.load()
    starttime = base_time.isel(times=0)
    # Count tracks within the specified period
    ntracks = len(starttime.where((starttime >= start_date) & (starttime <= end_date), drop=True))
    print(f'Number of tracks ({run_name}): {ntracks}')

    # Round the lat/lon to the nearest integer, 
    # This way the lat/lon are calculated at the precision of 1 degree (i.e., count on a 1x1 degree grid)
    # Note: counting using histogram2d only works on 1x1 degree grid (unique lat/lon are round to integer)
    rlat = ds['meanlat'].where((starttime >= start_date) & \
        (starttime <= end_date), drop=True).load().round().data
    rlon = ds['meanlon'].where((starttime >= start_date) & \
        (starttime <= end_date), drop=True).load().round().data

    rlat0 = ds['meanlat'].isel(times=0).where((starttime >= start_date) & \
        (starttime <= end_date), drop=True).load().round().data
    rlon0 = ds['meanlon'].isel(times=0).where((starttime >= start_date) & \
        (starttime <= end_date), drop=True).load().round().data

    # Total tracks count
    track_counts = count_unique_tracks(rlat, rlon, xbins, ybins)
    # # Start time
    # start_datetime = pd.to_datetime(sdate)

    # Calculate lat/lon bin center value
    xbins_c = xbins[:-1] + np.diff(xbins)/2.
    ybins_c = ybins[:-1] + np.diff(ybins)/2.

    # Write output file
    var_dict = {
        'track_counts': (['time', 'lat', 'lon'], np.expand_dims(track_counts, axis=0)),
    }
    coord_dict = {
        'time': (['time'], np.expand_dims(start_date, axis=0)),
        'lat': (['lat'], ybins_c),
        'lon': (['lon'], xbins_c),
        'lat_bnds': (['lat_bnds'], ybins),
        'lon_bnds': (['lon_bnds'], xbins),
    }
    gattr_dict = {
        'title': 'Track counts on grid',
        'start_date': sdate,
        'end_date': edate,
        'contact':'Zhe Feng, zhe.feng@pnnl.gov',
        'created_on':time.ctime(time.time()),
    }
    dsout = xr.Dataset(var_dict, coords=coord_dict, attrs=gattr_dict)
    dsout['lon'].attrs['long_name'] = 'Longitude grid center value'
    dsout['lon'].attrs['units'] = 'degree'
    dsout['lat'].attrs['long_name'] = 'Latitude grid center value'
    dsout['lat'].attrs['units'] = 'degree'
    dsout['lon_bnds'].attrs['long_name'] = 'Longitude grid bounds'
    dsout['lon_bnds'].attrs['units'] = 'degree'
    dsout['lat_bnds'].attrs['long_name'] = 'Latitude grid grid bounds'
    dsout['lat_bnds'].attrs['units'] = 'degree'
    dsout['track_counts'].attrs['long_name'] = 'Track counts on gridded map'
    dsout['track_counts'].attrs['units'] = 'count'

    fillvalue = np.nan
    # Set encoding/compression for all variables
    comp = dict(zlib=True, _FillValue=fillvalue, dtype='float32')
    encoding = {var: comp for var in dsout.data_vars}
    # Write output
    dsout.to_netcdf(path=output_filename, mode='w', format='NETCDF4', unlimited_dims='time', encoding=encoding)
    print(f'Output saved: {output_filename}')

    # import pdb; pdb.set_trace()

    # xbins5 = np.arange(0,361,5)
    # ybins5 = np.arange(-60,61,5)
    # ranges5 = [[min(ybins5), max(ybins5)], [min(xbins5), max(xbins5)]]
    # # ccs_initcount_obs, yedges, xedges = np.histogram2d(lat0_obs, lon0_obs, bins=[ybins5, xbins5], range=ranges5)
    # ccs_initcount_mod, yedges, xedges = np.histogram2d(lat0_mod, lon0_mod, bins=[ybins5, xbins5], range=ranges5)