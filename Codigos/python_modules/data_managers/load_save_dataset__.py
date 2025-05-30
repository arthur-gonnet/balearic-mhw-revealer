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


# USER INPUT - Bathy dataset path
bathy_GEBCO_path    = os.path.join(datos_location, "bathymetry", "Bathymetry_GEBCO_2023_IBERIAN.nc")
bathy_MFC_path      = os.path.join(datos_location, "bathymetry", "Bathymetry_MED-MFC_006_004_mask_bathy_BalearicIslands.nc")
bathy_MFC_med_path  = os.path.join(datos_location, "bathymetry", "Bathymetry_MED-MFC_006_004_mask_bathy.nc")

# USER INPUT - Dragonera insitu dataset path
dragonera_insitu_path           = os.path.join(datos_location, "puerto_del_estado", "dragonera_hist_temp.nc")
dragonera_insitu_marta_s_path   = os.path.join(datos_location, "puerto_del_estado", "dragonera_hist_temp_marta_s.nc")

# USER INPUT - CMEM-SST pathes
cmems_sst_folder_path   = os.path.join(datos_location, "CMEMS-SST", "MEDITERRANEAN", "SST-L4-REP-HR", "DATA-NEW", "DAILY")
cmems_sst_ds_pattern    = os.path.join(cmems_sst_folder_path, "{year}", "SST_MED_SST_L4_REP_OBSERVATIONS_010_021_y{year}m{month}.nc")

cmems_sst_trend_path    = os.path.join(datos_location, "CMEMS-SST", "processed", "cmems_sst_trend_1982_2023.nc")
cmems_sst_clim_pattern  = os.path.join(datos_location, "CMEMS-SST", "processed", "cmems_sst_clim_{region}_{year_start}_{year_end}.nc")

# USER INPUT - CMEM-MFC pathes
cmems_mfc_folder_path   = os.path.join(datos_location, "CMEMS-MFC", "MEDITERRANEAN", "REANALYSIS", "DATA", "DAILY", "BalearicIslands")
cmems_mfc_ds_pattern    = os.path.join(cmems_mfc_folder_path, "{year}", "TEMP_MEDSEA_MULTIYEAR_PHY_006_004_y{year}m{month}_BalearicIslands.nc")

# USER INPUT - MHWs datasets pathes
mhws_dataset_pattern    = os.path.join(datos_location, "mhws", "{dataset}_mhws_{region}_{clim_start}_{clim_end}.nc")


########################################################################################################################
##################################### CODES ############################################################################
########################################################################################################################

# # # # Bathymetry # # # #

def load_bathy(source: str="GEBCO") -> xr.Dataset:
    """
    Loads a bathymetry dataset using xarray (either GEBCO or MFC).

    Parameters
    ----------
    source: str, default="GEBCO"
        Source of the bathymetry dataset to load. Can be "GEBCO", "MFC" or "MFC_MED".
    
    Returns
    ----------
    ds_bathy: xr.Dataset
        The bathymetry dataset.
    """

    if source.lower() == "gebco":
        # Loading the dataset
        ds_bathy_GEBCO = load_nc_to_dataset(bathy_GEBCO_path, remove_bay_of_biscay=True)

        # Adjusting the depth parameter, name and value
        # The expected depth parameter should be positive at depth
        ds_bathy_GEBCO = ds_bathy_GEBCO.rename({
            "elevation": "depth"
        })
        ds_bathy_GEBCO = ds_bathy_GEBCO.where(ds_bathy_GEBCO['depth'] < 10)
        ds_bathy_GEBCO['depth'] = -ds_bathy_GEBCO['depth']

        # Printing the good news
        print("Loaded GEBCO bathymetry dataset.")

        # Returning the final dataset
        return ds_bathy_GEBCO
    
    elif source.lower() == "mfc":
        # Loading the dataset
        ds_bathy_MFC = load_nc_to_dataset(
            bathy_MFC_path,
            remove_bay_of_biscay=False,
            drop_variables=["deptho_lev", "mask"]
        )

        # Adjusting the parameters, name and value
        ds_bathy_MFC = ds_bathy_MFC.drop_dims('depth')
        ds_bathy_MFC = ds_bathy_MFC.rename({
            "longitude": "lon",
            "latitude": "lat",
            "deptho": "depth"
        })


        # Printing the good news
        print("Loaded MFC bathymetry dataset.")

        # Returning the final dataset
        return ds_bathy_MFC
    
    elif source.lower() == "mfc_med":
        # Loading the dataset
        ds_bathy_MFC = load_nc_to_dataset(bathy_MFC_med_path, remove_bay_of_biscay=True)

        # Printing the good news
        print("Loaded all Mediterranean MFC bathymetry dataset.")

        # Returning the final dataset
        return ds_bathy_MFC
    
    else:
        print(f"Bathy dataset not recognized ({source}), please choose 'gebco' or 'mfc'.")

# # # # Insitu # # # #

def load_dragonera_insitu(
        time_selector: str | slice | None = None,

        nightime_data: bool = False,
        daily_avg: bool = False,
        drop_depth: bool = False,

        load_marta_s: bool = False,
) -> xr.Dataset:
    """
    Loads the Dragonera insitu dataset using xarray with optional temporal subsetting.

    Parameters
    ----------
    time_selector: str | slice[str] | None, optional
        Time selection applied using xarray's `.sel()`. Can be a string (e.g., "1993-01-21") or a slice
        (e.g., slice("1993-01-21", "1993-01-25")). If None, no time filtering is applied.
    
    nightime_data: bool, default=False
        If True, only the data between 9 p.m. and 6 a.m. (excluded) will be kept in the dataset.
        This is meant to match the CMEMS-SST time criterion.
    
    daily_avg: bool, default=False
        If True, daily averages the dataset. If nan value is present into a day, the resulting average of that day will be nan.
    
    drop_depth: bool, default=False
        If True, the depth coordinate, set as 3m, will be discarded.
    
    load_marta_s: bool, default=False
        If True, the dataset of Marta is loaded instead of the PdE dataset.
    
    Returns
    ----------
    ds_dragonera_insitu: xr.Dataset
        The Dragonera insitu dataset.
    """

    # Choosing either loading marta's file or normal file
    filepath = dragonera_insitu_marta_s_path if load_marta_s else dragonera_insitu_path
    
    # Loading the dataset
    ds_dragonera_insitu = xr.open_dataset(filepath)

    # Applying time selection
    if not time_selector is None:
        ds_dragonera_insitu = ds_dragonera_insitu.sel(time=time_selector, method=(None if type(time_selector) == slice else 'nearest'))

    # Adding coordinates to dataset
    if drop_depth:
        ds_dragonera_insitu = ds_dragonera_insitu.assign_coords({"lon": 2.10, "lat": 39.56})
    else:
        ds_dragonera_insitu = ds_dragonera_insitu.assign_coords({"lon": 2.10, "lat": 39.56, "depth": 3})

    # Keep only the nighttime data
    if nightime_data:
        ds_dragonera_insitu = ds_dragonera_insitu.isel(time=(
            (ds_dragonera_insitu.time.dt.hour < 6) |
            (ds_dragonera_insitu.time.dt.hour >= 21)
        ))

    # Resampling if necessary
    if daily_avg:
        # TODO : Should load precalculated dataset

        # Averages values
        ds_dragonera_insitu = ds_dragonera_insitu.resample(time="1D").mean(skipna=False)
        # # Either do the daily mean centered at midnight or at midday
        # if daily_avg_center_0am:
        #     ds_dragonera_insitu = ds_dragonera_insitu.resample(time="1D", offset="12h").mean()
        
        #     # Moving time so that value is on the exact day still
        #     offset = pd.tseries.frequencies.to_offset("12h")
        #     ds_dragonera_insitu["time"] = ds_dragonera_insitu.get_index("time") + offset
        
        # else:

    # Saving unit as attribute
    ds_dragonera_insitu.time.attrs["reference"] = "GMT"
    ds_dragonera_insitu.T.attrs["unit"] = "°C"

    # Printing the good news
    details = []
    if nightime_data: details.append("nightime")
    if daily_avg: details.append("daily averaged")
    details = (" (" + ', '.join(details) + ")") if details != [] else ''
    
    print(f"Loaded Dragonera insitu dataset{details}.")

    # Returning the final dataset
    return ds_dragonera_insitu

# # # # CMEMS-SST # # # #

def load_cmems_sst(
        years: list[int | str] = range(1982, 2024),
        months: list[int | str] = range(1, 13),

        time_selector: str | slice | None = None,
        lon_selector: float | slice | None = None,
        lat_selector: float | slice | None = None,
        region_selector: str | None = None,
        
        chunks: int | dict | str | None = 'auto',

        only_sst: bool = True,
) -> xr.Dataset:
    """
    Loads the CMEMS-SST dataset using xarray. By default, the entire dataset is loaded, but optional spatio-temporal subsettings are available
    for faster loading times.

    Parameters
    ----------
    years: list[int | str], default=['*'], optional
        Years to load from the dataset, as integers or strings (e.g., range(1983, 1985) or ["1983"]).
        Use ['*'] to load all available years (1982-2023) as it uses a glob pattern.
        This parameter defines the files to load, so loading time is proportional with this range size.
    
    months: list[int | str], default=['*'], optional
        Months to load from the dataset, as integers or strings (e.g., range(1, 5) or ["1"]).
        Use ['*'] to load all months (1-12) as it uses a glob pattern.
        This parameter defines the files to load, so loading time is proportional with this range size.
    
    time_selector: str | slice[str] | None, optional
        Time selection applied using xarray's `.sel()`. Can be a string (e.g., "1993-01-21") or a slice
        (e.g., slice("1993-01-21", "1993-01-25")). If None, no time filtering is applied.
    
    lon_selector: float | slice[float] | None, optional
        Longitude selector applied using xarray's `.sel()`. Accepts a float or a slice.
    
    lat_selector: float | slice[float] | None, optional
        Latitude selector applied using xarray's `.sel()`. Accepts a float or a slice.
    
    chunks: int | dict | 'auto' | None, default='auto'
        xarray parameter.
        Dictionary with keys given by dimension names and values given by chunk sizes.
        In general, these should divide the dimensions of each dataset. If int, chunk
        each dimension by ``chunks``. By default, chunks will be chosen to load entire
        input files into memory at once. This has a major impact on performance.
    
    only_sst: bool, default=False
        If True, all other variables than SST will be discarded.
    
    Returns
    ----------
    ds_cmems_sst: xr.Dataset
        The CMEMS-SST dataset.
    """

    # Selecting the files to load
    files: list[str] = []

    for year in years:
        if year != '*' and (int(year) < 1982 or int(year) > 2023):
            print(f"Incorrect year selection for CMEMS-SST dataset, got {years}, expected a year range between 1982 and 2023")
        
        for month in months:
            if month != '*' and (int(month) < 1 or int(month) > 12):
                print(f"Incorrect month selection for CMEMS-SST dataset, got {months}, expected a month range between 1 and 12")
            
            month_str: str = str(month) if (month == '*' or month > 9) else ('0'+str(month))
            pattern: str = cmems_sst_ds_pattern.format(year=year, month=month_str)
            files.extend(glob.glob(pattern))
    
    if region_selector == "balears":
        lon_selector = slice(-0.9, 5.1)
        lat_selector = slice(37.6, 41.1)

    # Preprocess applied before concatenating the datasets
    def preprocess(ds: xr.Dataset) -> xr.Dataset:
        """
        The preprocess applies dataset operations before concatenating the multiple datasets.
        Here, performing spatial selection.
        """
        
        # Remove bay of Biscay
        ds = ds.where(((ds.lon > 0) | (ds.lat < 42)), drop=True)

        # Spatial selection
        if not lon_selector is None:
            # Depending if lon_selector is slice or float, use different method
            ds = ds.sel(lon=lon_selector, method=(None if type(lon_selector) == slice else 'nearest'))
        
        if not lat_selector is None:
            ds = ds.sel(lat=lat_selector, method=(None if type(lat_selector) == slice else 'nearest'))
        
        return ds

    # Options for loading the dataset
    drop_vars = ["analysis_error", "mask", "sea_ice_fraction"] if only_sst else None


    # LOADING THE DATASET !
    ds_cmems_sst = xr.open_mfdataset(
        paths = files,
        preprocess = preprocess,
        drop_variables = drop_vars,
        decode_times = True,
        chunks = chunks,
    )


    # Time selection (some issues occur when selecting in preprocess)
    if not time_selector is None:
        ds_cmems_sst = ds_cmems_sst.sel(time=time_selector)
    
    # Uniformying variables name to T for temperature
    ds_cmems_sst = ds_cmems_sst.rename({"analysed_sst": "T"})

    # Changing unit from K to °C
    ds_cmems_sst["T"] = ds_cmems_sst.T - 273.15
    ds_cmems_sst.T.attrs["unit"] = "°C"

    # Printing the good news
    print(f"Loaded CMEMS-SST dataset.")
    
    # Returning the final dataset
    return ds_cmems_sst

def save_cmems_sst_trend(
        ds_cmems_sst_trend: xr.Dataset,
        progress_bar: bool = True,
        profilers: bool = False
):
    """
    Saves a CMEMS-SST trend dataset using xarray.

    Parameters
    ----------
    ds_cmems_sst_trend: xr.Dataset, 
        CMEMS-SST trend dataset to save.
    
    progress_bar: bool, default=True
        If True, shows a progress bar when saving the dataset.
    
    profilers: bool, default=True
        If True, shows the profilers after that the dataset is saved.
    """
    
    # Printing the good news
    print(f"Saving CMEMS-SST trend dataset to {cmems_sst_trend_path}")

    # Adding metadata
    ds_cmems_sst_trend.sst_trend.attrs["unit"] = "°C/year"
    ds_cmems_sst_trend.attrs["title"] = "Linear SST trend over the period 1982-2023 in the Balearic sea using " \
                                "Mediterranean Sea SST Analysis (cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021)"
    ds_cmems_sst_trend.attrs["acknowledgment"] = "Generated/provided by Copernicus Marine Service and CNR - ISMAR ROME."

    save_dataset_to_nc(
        ds = ds_cmems_sst_trend,
        file_path = cmems_sst_trend_path,
        progress_bar = progress_bar,
        profilers = profilers
    )

    print(f" -> Saved!")

def load_cmems_sst_trend() -> xr.Dataset:
    """
    Loads the CMEMS-SST trend dataset using xarray.
    
    Returns
    ----------
    ds_cmems_sst_trend: xr.Dataset
        The loaded CMEMS-SST trend dataset.
    """

    # Loading the dataset
    ds_cmems_sst_trend = load_nc_to_dataset(cmems_sst_trend_path, remove_bay_of_biscay=True)

    # Printing the good news
    print("Loaded CMEMS-SST trend.")
    
    # Returning the final dataset
    return ds_cmems_sst_trend

def save_cmems_sst_clim(
        ds_cmems_sst_clim: xr.Dataset,

        clim_year_start: int = 1987,
        clim_year_end: int = 2021,
        region: str = "balears",

        progress_bar: bool = True,
        profilers: bool = False,
):
    """
    Saves a CMEMS-SST climatology dataset using xarray.

    Parameters
    ----------
    ds_cmems_sst_clim: xr.Dataset, 
        CMEMS-SST climatology dataset to save.

    clim_year_start: int, default=1987
        Year that the climatology period started. Will influence the output file name.

    clim_year_end: int, default=2021
        Year that the climatology period ended. Will influence the output file name.
    
    progress_bar: bool, default=True
        If True, shows a progress bar when saving the dataset.
    
    profilers: bool, default=True
        If True, shows the profilers after that the dataset is saved.
    """

    # Choosing file path corresponding to climatology period
    cmems_sst_clim_path = cmems_sst_clim_pattern.format(region=region, year_start=clim_year_start, year_end=clim_year_end)
    
    # Printing the good news
    print(f"Saving CMEMS-SST climatology dataset to {cmems_sst_clim_path}")

    # Adding metadata
    ds_cmems_sst_clim.attrs["title"] = "Climatology over the period {clim_year_start}-{clim_year_end} in the Balearic sea using " \
                                "Mediterranean Sea SST Analysis (cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021)"
    ds_cmems_sst_clim.attrs["acknowledgment"] = "Generated/provided by Copernicus Marine Service and CNR - ISMAR ROME."

    save_dataset_to_nc(
        ds = ds_cmems_sst_clim,
        file_path = cmems_sst_clim_path,
        progress_bar = progress_bar,
        profilers = profilers
    )

    print(f" -> Saved!")

def load_cmems_sst_clim(
        clim_year_start: int=1987,
        clim_year_end: int=2021,
        region: str = "balears",
) -> xr.Dataset:
    """
    Loads the CMEMS-SST climatology dataset using xarray.

    Parameters
    ----------
    clim_year_start: int, default=1987
        Desired year of climatology period start. Defines the input file name.

    clim_year_end: int, default=2021
        Desired year of climatology period end. Defines the input file name.
    
    Returns
    ----------
    ds_cmems_sst_clim: xr.Dataset
        The loaded CMEMS-SST climatology dataset.
    """

    # Choosing file path corresponding to climatology period
    cmems_sst_clim_path = cmems_sst_clim_pattern.format(region=region, year_start=clim_year_start, year_end=clim_year_end)

    # Loading the dataset
    ds_cmems_sst_clim = load_nc_to_dataset(
        file_path = cmems_sst_clim_path,
        remove_bay_of_biscay = True
    )

    # Printing the good news
    print("Loaded CMEMS-SST climatology.")
    
    # Returning the final dataset
    return ds_cmems_sst_clim

# # # # CMEMS-MFC # # # #

def load_cmems_mfc(
        years: list[int|str] = range(1987, 2023),
        months: list[int|str] = range(1, 13),

        time_selector: str|slice|None = None,
        lon_selector: float|slice|None = None,
        lat_selector: float|slice|None = None,
        depth_selector: float|slice|None = None,
        region_selector: str | None = None,

        move_to_0am: bool = True,

        drop_vars = None,
        
        chunks: int | dict | str | None = None,
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


    # LOADING THE DATASET !
    ds_cmems_mfc = xr.open_mfdataset(
        files,
        preprocess=preprocess,
        drop_variables=drop_vars,
        decode_times=True,
        chunks=chunks,
        # combine="nested"
    )

    # Applying time selection (some issues occur when selecting in preprocess)
    if not time_selector is None:
        ds_cmems_mfc = ds_cmems_mfc.sel(time=time_selector)
    
    # Using same variable name across everydataset
    if not "thetao" in drop_vars:
        ds_cmems_mfc = ds_cmems_mfc.rename({"thetao": "T"})
    
    # Resampling if necessary
    if move_to_0am:
        ds_cmems_mfc["time"] = ds_cmems_mfc.time.dt.floor("1D")

    # Adding custom attributes to dataset
    if not "thetao" in drop_vars:
        ds_cmems_mfc.T.attrs["unit"] = "°C"

    # Printing the good news
    print("Loaded CMEMS-MFC dataset.")

    # Returning the final dataset
    return ds_cmems_mfc

# # # # Generic # # # #

def save_dataset_to_nc(
        ds: xr.Dataset,
        file_path: str,
        progress_bar: bool = True,
        profilers: bool = False,
) -> xr.Dataset:
    if progress_bar:
        with ProgressBar():
            ds = ds.compute()
            ds.to_netcdf(file_path)
    
    elif progress_bar and profilers:
        with ProgressBar(), ResourceProfiler() as rprof, Profiler() as prof:
            ds = ds.compute()
            ds.to_netcdf(file_path)
        
        from datetime import datetime
        path = os.path.join(codigos_location, ".dask_profiles", datetime.now().isoformat() + ".html")
        visualize([prof, rprof], path, show=True, save=True)

    else:
        ds = ds.compute()
        ds.to_netcdf(file_path)

    return ds

def load_nc_to_dataset(
        file_path,
        remove_bay_of_biscay: bool = False,
        **kwargs
):
    # Loading the dataset
    ds = xr.open_dataset(
        file_path,
        decode_times=True,
        **kwargs
    )

    # Remove bay of Biscay
    if remove_bay_of_biscay:
        ds = ds.where(((ds.lon > 0) | (ds.lat < 42)), drop=True)
    
    # Returning the final dataset
    return ds


def save_mhws_dataset(
        ds_mhws: xr.Dataset,

        original_dataset: str,
        clim_period: list[int],
        region: str = "balears",

        progress_bar: bool = True,
        profilers: bool = False
):
    """
    Saves a MHWs dataset using xarray.

    Parameters
    ----------
    ds_mhws: xr.Dataset
        MHWs dataset to save.

    original_dataset: str 
        Dataset from which the MHWs computations were performed. Can be "cmems_sst" or "mfc_50m".

    region: str, default="balears"
        Region in which the MHWs computations were performed. Can be "balears" or "med"

    clim_period: list[int]
        Climatology period used for the MHWs computations.
    
    progress_bar: bool, default=True
        If True, shows a progress bar when saving the dataset.
    
    profilers: bool, default=True
        If True, shows the profilers after that the dataset is saved.
    """
    
    mhws_dataset_path = mhws_dataset_pattern.format(dataset=original_dataset, region=region, clim_start=clim_period[0], clim_end=clim_period[1])

    # Printing the good news
    print(f"Saving MHWs dataset to {mhws_dataset_path}")

    ds_mhws = save_dataset_to_nc(
        ds = ds_mhws,
        file_path = mhws_dataset_path,
        progress_bar = progress_bar,
        profilers = profilers
    )

    print(f" -> Saved!")

    return ds_mhws

def load_mhws_dataset(
        original_dataset: str,
        clim_period: tuple[int, int],
        region: str = "balears",
):
    """
    Loads the MHWs dataset using xarray.

    Parameters
    ----------
    original_dataset: str 
        Dataset from which the MHWs computations were performed. Can be "cmems_sst" or "mfc_50m".

    region: str, default="balears"
        Region in which the MHWs computations were performed. Can be "balears" or "med"

    clim_period: list[int]
        Climatology period used for the MHWs computations.
    
    Returns
    ----------
    ds_mhws: xr.Dataset
        The loaded MHWs dataset.
    """

    # Choosing file path corresponding to climatology period
    mhws_dataset_path = mhws_dataset_pattern.format(dataset=original_dataset,region=region, clim_start=clim_period[0], clim_end=clim_period[1])

    # Loading the dataset
    ds_mhws = load_nc_to_dataset(
        file_path = mhws_dataset_path,
        remove_bay_of_biscay = True
    )

    # Printing the good news
    print("Loaded MHWs dataset.")
    
    # Returning the final dataset
    return ds_mhws
