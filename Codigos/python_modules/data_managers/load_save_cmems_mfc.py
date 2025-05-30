########################################################################################################################
########################################################################################################################

###### USER NOTES:
# This script gather all the codes to load and save datasets using xarray. 

########################################################################################################################
##################################### IMPORTS ##########################################################################
########################################################################################################################

import os
import glob as glob

import pandas as pd
import xarray as xr

from dask.diagnostics import ProgressBar, ResourceProfiler, Profiler, visualize

codigos_location = os.path.abspath(os.getcwd())
datos_location = os.path.join(codigos_location, "..", "Datos")

########################################################################################################################
##################################### USER INPUT #######################################################################
########################################################################################################################

# USER INPUT - CMEM-MFC pathes
cmems_mfc_folder_path   = os.path.join(datos_location, "CMEMS-MFC", "MEDITERRANEAN", "REANALYSIS", "DATA", "DAILY", "BalearicIslands")
cmems_mfc_ds_pattern    = os.path.join(cmems_mfc_folder_path, "{year}", "TEMP_MEDSEA_MULTIYEAR_PHY_006_004_y{year}m{month}_BalearicIslands.nc")


########################################################################################################################
##################################### CODES ############################################################################
########################################################################################################################

def load_cmems_mfc(
        years: list[int|str] = range(1987, 2023),
        months: list[int|str] = range(1, 13),

        time_selector: str|slice|None = None,
        lon_selector: float|slice|None = None,
        lat_selector: float|slice|None = None,
        depth_selector: float|slice|None = None,
        region_selector: str | None = None,

        move_to_0am: bool = True,
        only_botT: bool = False,

        drop_vars = [],
        
        chunks: int | dict | str | None = "auto",
) -> xr.Dataset:
    """
    Loads the CMEMS-MFC dataset using xarray. By default, the entire dataset is loaded, but optional spatio-temporal subsettings are available
    for faster loading times.

    Parameters
    ----------
    years: list[int|str], default=['*']
        Years to load from the dataset, as integers or strings (e.g., range(1993, 1995) or ["1993"]).
        Use ['*'] to load all available years (1987-2022) as it uses a glob pattern.
    
    months: list[int|str], default=['*']
        Months to load from the dataset, as integers or strings (e.g., range(1, 5) or ["1"]).
        Use ['*'] to load all months (1-12) as it uses a glob pattern.
    
    time_selector: str | slice[str] | None, default=None, optional
        Time selection applied using xarray's `.sel()`. Can be a string (e.g., "1993-01-21") or a slice
        (e.g., slice("1993-01-21", "1993-01-25")). If None, no time filtering is applied.
    
    lon_selector: float | slice[float] | None, default=None, optional
        Longitude selector applied using xarray's `.sel()`. Accepts a float or a slice. If None, no longitude filtering is applied.
    
    lat_selector: float | slice[float] | None, default=None, optional
        Latitude selector applied using xarray's `.sel()`. Accepts a float or a slice. If None, no latitude filtering is applied.
    
    chunks: int | dict | 'auto' | None, optional
        xarray parameter.
        Dictionary with keys given by dimension names and values given by chunk sizes.
        In general, these should divide the dimensions of each dataset. If int, chunk
        each dimension by ``chunks``. By default, chunks will be chosen to load entire
        input files into memory at once. This has a major impact on performance:.
    
    Returns
    ----------
    ds_cmems_mfc: xr.Dataset
        The loaded CMEMS-MFC dataset with optional spatio-temporal subsetting.
    """


    # Selecting the files to load
    files: list[str] = []

    for year in years:
        if year != '*' and (int(year) < 1982 or int(year) > 2023):
            print(f"Incorrect year selection for CMEMS-MFC dataset, got {years}, expected a year range between 1982 and 2023")
        
        for month in months:
            if month != '*' and (int(month) < 1 or int(month) > 12):
                print(f"Incorrect month selection for CMEMS-MFC dataset, got {months}, expected a month range between 1 and 12")
            
            month_str: str = str(month) if month == '*' or month > 9 else '0'+str(month)
            pattern: str = cmems_mfc_ds_pattern.format(year=year, month=month_str)
            files.extend(glob.glob(pattern))
    
    if region_selector == "balears":
        lon_selector = slice(-0.9, 5.1)
        lat_selector = slice(37.6, 41.1)
    
    # Preprocess applied before concatenating the datasets
    def preprocess(ds: xr.Dataset) -> xr.Dataset:
        """ This preprocess selects the data while loading, which makes an efficient loading """
        
        # Remove bay of Biscay
        if not region_selector == "balears":
            ds = ds.where(((ds.lon > 0) | (ds.lat < 42)), drop=True)
        
        # Spatial selection
        if not lon_selector is None:
            # Depending if lon_selector is slice or float, use different method
            ds = ds.sel(lon=lon_selector, method=(None if type(lon_selector) == slice else 'nearest'))
        
        if not lat_selector is None:
            ds = ds.sel(lat=lat_selector, method=(None if type(lat_selector) == slice else 'nearest'))
        
        if not depth_selector is None:
            ds = ds.sel(depth=depth_selector, method=(None if type(depth_selector) == slice else 'nearest'))
        
        return ds

    # Acquiring only bottom temperature, those don't matter
    if only_botT:
        drop_vars = ["thetao", "depth"]


    # LOADING THE DATASET !
    ds_cmems_mfc = xr.open_mfdataset(
        files,
        preprocess=preprocess,
        drop_variables=drop_vars,
        decode_times=True,
        chunks=chunks,
    )

    # Applying time selection (some issues occur when selecting in preprocess)
    if not time_selector is None:
        ds_cmems_mfc = ds_cmems_mfc.sel(time=time_selector)
    
    # Resampling if necessary
    if move_to_0am:
        ds_cmems_mfc["time"] = ds_cmems_mfc.time.dt.floor("1D")

    if only_botT:
        # Compute the depth
        ...
    
    # Using same variable name across everydataset
    if not "thetao" in drop_vars:
        ds_cmems_mfc = ds_cmems_mfc.rename({"thetao": "T"})

    # Adding custom attributes to dataset
    if not "thetao" in drop_vars:
        ds_cmems_mfc.T.attrs["unit"] = "°C"

    # Printing the good news
    print("Loaded CMEMS-MFC dataset.")

    # Returning the final dataset
    return ds_cmems_mfc
