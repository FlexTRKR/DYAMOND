import numpy as np
import glob, os, sys
import xarray as xr
import pandas as pd
import time
import yaml

#--------------------------------------------------------------------------------------------
def calc_pdf(datafiles, outfile, lonbox, latbox, landmask_file):
    status = 0

    # Read data
    print(f'Reading input data ...')
    dropvars_list = ['numclouds','pcptracknumber']
    ds = xr.open_mfdataset(datafiles, concat_dim='time', combine='nested', drop_variables=dropvars_list)
    print(f'Finish reading data.')
    # Change lon from [0,360] to [-180,180]
    # latitude = ds.latitude.isel(time=0).load()
    # longitude = ((ds.longitude.isel(time=0).load() - 180) % 360) - 180
    # ds.coords['lon'] = ((ds.lon.load() - 180) % 360) - 180
    # ds['longitude'] = longitude
    # ds['latitude'] = latitude
    # lon = ds.longitude.isel(time=0).load()
    # lat = ds.latitude.isel(time=0).load()

    # Read landmask
    dslm = xr.open_dataset(landmask_file)
    landmask = dslm.landseamask
    # Replace landmask lon/lat coordinates to make sure they are the same (ignore round off error)
    landmask['lon'] = ds.lon.values
    landmask['lat'] = ds.lat.values

    longitude = ds.longitude.isel(time=0).load()
    latitude = ds.latitude.isel(time=0).load()

    # Create a mask for the region
    land_mask = (longitude >= lonbox[0]) & (longitude <= lonbox[1]) & (latitude >= latbox[0]) & (latitude <= latbox[1]) & (landmask <= 10)
    ocean_mask = (longitude >= lonbox[0]) & (longitude <= lonbox[1]) & (latitude >= latbox[0]) & (latitude <= latbox[1]) & (landmask >= 95)
    # Add mask as coordinates to the dataset
    ds.coords['landmask'] = (('lat', 'lon'), land_mask.data)
    ds.coords['oceanmask'] = (('lat', 'lon'), ocean_mask.data)

    # Set up the rain rate bins
    rrbins = np.arange(1, 201, 1)
    # rrbins = np.logspace(np.log10(0.01), np.log10(100.0), 100)
    rr_range = (np.min(rrbins), np.max(rrbins))

    # Select monsoon region over land
    totpcp_land = ds['precipitation'].where(ds.landmask == 1, drop=True)
    totpcp_land_pdf, bins = np.histogram(totpcp_land, bins=rrbins, range=rr_range, density=False)
    del totpcp_land

    # For MCS precipitation, use cloudtracknumber to mask (replace non-MCS area with 0 for averaging purpose), then discard the region outside
    mcspcp_land = ds['precipitation'].where(ds['cloudtracknumber'] > 0).where(ds.landmask == 1, drop=True)
    mcspcp_land_pdf, bins = np.histogram(mcspcp_land, bins=rrbins, range=rr_range, density=False)
    del mcspcp_land

    # Non-MCS deep convection (cloudnumber > 0: CCS & cloudtracknumber == NaN: non-MCS)
    idcpcp_land = ds['precipitation'].where((ds['cloudnumber'] > 0) & (np.isnan(ds['cloudtracknumber']))).where(ds.landmask == 1, drop=True)
    idcpcp_land_pdf, bins = np.histogram(idcpcp_land, bins=rrbins, range=rr_range, density=False)
    del idcpcp_land

    # Congestus (Tb between CCS (241K) and 310 K, rain rate > 0.5 mm/h)
    congpcp_land = ds['precipitation'].where(np.isnan(ds['cloudnumber']) & (ds.tb < 310) & (ds['precipitation'] > 0.5)).where(ds.landmask == 1, drop=True)
    # shcupcp_land = ds['precipitation'].where(np.isnan(ds['cloudnumber'])).where(ds.landmask == 1, drop=True)
    congpcp_land_pdf, bins = np.histogram(congpcp_land, bins=rrbins, range=rr_range, density=False)
    del congpcp_land

    # Select monsoon region over ocean
    totpcp_ocean = ds['precipitation'].where(ds.oceanmask == 1, drop=True)
    totpcp_ocean_pdf, bins = np.histogram(totpcp_ocean, bins=rrbins, range=rr_range, density=False)
    del totpcp_ocean

    # For MCS precipitation, use cloudtracknumber to mask (replace non-MCS area with 0 for averaging purpose), then discard the region outside
    mcspcp_ocean = ds['precipitation'].where(ds['cloudtracknumber'] > 0).where(ds.oceanmask == 1, drop=True)
    mcspcp_ocean_pdf, bins = np.histogram(mcspcp_ocean, bins=rrbins, range=rr_range, density=False)
    del mcspcp_ocean

    # Non-MCS deep convection (cloudnumber > 0: CCS & cloudtracknumber == NaN: non-MCS)
    idcpcp_ocean = ds['precipitation'].where((ds['cloudnumber'] > 0) & (np.isnan(ds['cloudtracknumber']))).where(ds.oceanmask == 1, drop=True)
    idcpcp_ocean_pdf, bins = np.histogram(idcpcp_ocean, bins=rrbins, range=rr_range, density=False)
    del idcpcp_ocean

    # Congestus (Tb between CCS (241K) and 310 K, rain rate > 0.5 mm/h)
    congpcp_ocean = ds['precipitation'].where(np.isnan(ds['cloudnumber']) & (ds.tb < 310) & (ds['precipitation'] > 0.5)).where(ds.oceanmask == 1, drop=True)
    congpcp_ocean_pdf, bins = np.histogram(congpcp_ocean, bins=rrbins, range=rr_range, density=False)
    del congpcp_ocean  


    # Define xarray output dataset
    print('Writing output to netCDF file ...')
    var_dict = {
        'total_land': (['bins'], totpcp_land_pdf),
        'mcs_land': (['bins'], mcspcp_land_pdf),
        'idc_land': (['bins'], idcpcp_land_pdf),
        'congestus_land': (['bins'], congpcp_land_pdf),
        'total_ocean': (['bins'], totpcp_ocean_pdf),
        'mcs_ocean': (['bins'], mcspcp_ocean_pdf),
        'idc_ocean': (['bins'], idcpcp_ocean_pdf),
        'congestus_ocean': (['bins'], congpcp_ocean_pdf),
    }
    coord_dict = {'bins': (['bins'], rrbins[:-1])}
    gattr_dict = {
        'title': 'Precipitation PDF by types',
        'lonbox':lonbox,
        'latbox':latbox,
        'contact':'Zhe Feng, zhe.feng@pnnl.gov',
        'created_on':time.ctime(time.time()),
    }
    dsout = xr.Dataset(var_dict, coords=coord_dict, attrs=gattr_dict)

    dsout.bins.attrs['long_name'] = 'Rain rate bins'
    dsout.bins.attrs['units'] = 'mm/h'
    dsout.total_land.attrs['long_name'] = 'Land total precipitation'
    dsout.total_land.attrs['units'] = 'count'
    dsout.mcs_land.attrs['long_name'] = 'Land MCS precipitation'
    dsout.mcs_land.attrs['units'] = 'count'
    dsout.idc_land.attrs['long_name'] = 'Land isolated deep convection precipitation'
    dsout.idc_land.attrs['units'] = 'count'
    dsout.congestus_land.attrs['long_name'] = 'Land congestus precipitation'
    dsout.congestus_land.attrs['units'] = 'count'
    dsout.total_ocean.attrs['long_name'] = 'Ocean total precipitation'
    dsout.total_ocean.attrs['units'] = 'count'
    dsout.mcs_ocean.attrs['long_name'] = 'Ocean MCS precipitation'
    dsout.mcs_ocean.attrs['units'] = 'count'
    dsout.idc_ocean.attrs['long_name'] = 'Ocean isolated deep convection precipitation'
    dsout.idc_ocean.attrs['units'] = 'count'
    dsout.congestus_ocean.attrs['long_name'] = 'Ocean congestus precipitation'
    dsout.congestus_ocean.attrs['units'] = 'count'

    fillvalue = np.nan
    # Set encoding/compression for all variables
    comp = dict(zlib=True, dtype='float')
    encoding = {var: comp for var in dsout.data_vars}
    # Write to file
    dsout.to_netcdf(path=outfile, mode='w', format='NETCDF4', encoding=encoding)
    print('Output saved as: ', outfile)

    status = 1

    return status



if __name__ == "__main__":
    
    # Get runtime inputs
    config_file = sys.argv[1]
    run_name = sys.argv[2]

    # Get inputs from configuration file
    stream = open(config_file, 'r')
    config = yaml.full_load(stream)
    
    start_datetime = config['start_datetime']
    end_datetime = config['end_datetime']
    region = config['region']
    lon_bounds = config['lon_bounds']
    lat_bounds = config['lat_bounds']
    indir = config['indir']
    outdir = config['outdir']
    landmask_file = config['landmask_file']

    indir = f'{indir}{run_name}/mcstracking/20200120.0000_20200301.0000/'
    outdir = f'{outdir}{run_name}/stats/'

    # Generate time marks within the start/end datetime
    file_datetimes = pd.date_range(start=start_datetime, end=end_datetime, freq='1D').strftime('%Y%m%d')
    # Find all files from these dates
    datafiles = []
    for tt in range(0, len(file_datetimes)):
        datafiles.extend(sorted(glob.glob(f'{indir}mcstrack_{file_datetimes[tt]}*.nc')))
    print(f'Number of files: {len(datafiles)}')

    # Output filename
    outfile = f'{outdir}mcs_rainrate_hist_{file_datetimes[0]}_{file_datetimes[-1]}_{region}.nc'

    # Call function
    status = calc_pdf(datafiles, outfile, lon_bounds, lat_bounds, landmask_file)
