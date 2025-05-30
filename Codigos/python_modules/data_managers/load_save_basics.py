########################################################################################################################
########################################################################################################################

###### USER NOTES:
# This script gather all the codes to load and save datasets using xarray. 

## Functions description
#
# - load_bathy(...):
#       Loads a bathymetry dataset using xarray (either GEBCO or MFC).
#
# - load_dragonera_insitu(...):
#       Loads the Dragonera insitu dataset using xarray.
#
# - save_dataset_to_nc(...):
#       Basic function to save a dataset using xarray.
#
# - load_nc_to_dataset(...):
#       Basic function to load a dataset using xarray.

########################################################################################################################
##################################### IMPORTS ##########################################################################
########################################################################################################################

# Basic imports
import os
import glob as glob

# Advanced imports
import xarray as xr

from dask.diagnostics import ProgressBar, ResourceProfiler, Profiler, visualize

# Load relative file paths
# These locations must be accurate to load the data
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


########################################################################################################################
##################################### CODES ############################################################################
########################################################################################################################


def load_bathy(source: str="GEBCO", drop_vars="all", region_selector:str | None=None) -> xr.Dataset:
    """
    Loads a bathymetry dataset using xarray (either GEBCO or MFC).

    Parameters
    ----------
    source: str, default="GEBCO"
        Source of the bathymetry dataset to load. Can be `"GEBCO"`, `"MFC"` or `"MFC_MED"`.
    
    Returns
    ----------
    ds_bathy: xr.Dataset
        The loaded bathymetry dataset.
    """

    if source.lower() == "gebco":
        # Loading the dataset
        ds_bathy_GEBCO = load_nc_to_dataset(
            bathy_GEBCO_path,
            name_dict={
                "elevation": "depth"
            },
            remove_bay_of_biscay=True,
            region_selector=region_selector
        )

        # Adjusting the depth parameter, name and value
        # The expected depth parameter should be positive at depth
        ds_bathy_GEBCO = ds_bathy_GEBCO.where(ds_bathy_GEBCO['depth'] < 10)
        ds_bathy_GEBCO['depth'] = -ds_bathy_GEBCO['depth']

        # Printing the good news
        print("Loaded GEBCO bathymetry dataset.")

        # Returning the final dataset
        return ds_bathy_GEBCO
    
    elif source.lower() == "mfc" or source.lower() == "mfc_med":
        # Loading the dataset
        ds_bathy_MFC = load_nc_to_dataset(
            bathy_MFC_path if source.lower() == "mfc" else bathy_MFC_med_path,
            name_dict = {
                "longitude": "lon",
                "latitude": "lat",
            },
            remove_bay_of_biscay=False,
            region_selector=region_selector,
            drop_variables=["deptho_lev", "mask"] if drop_vars=="all" else drop_vars
        )

        # Adjusting the parameters, name and value
        ds_bathy_MFC = ds_bathy_MFC.drop_dims('depth')
        ds_bathy_MFC = ds_bathy_MFC.rename({
            "deptho": "depth"
        })

        # Printing the good news
        print(f"Loaded MFC bathymetry dataset{' (all med)' if source.lower() == 'mfc_med' else ''}.")

        # Returning the final dataset
        return ds_bathy_MFC
    
    else:
        print(f"Bathymetry dataset not recognized ({source}), please choose \"GEBCO\", \"MFC\" or \"MFC_MED\".")

def load_dragonera_insitu(
        # Xarray selectors
        time_selector: str | slice | None = None,

        # Dataset operations
        nightime_data: bool = False,
        daily_avg: bool = False,
        drop_depth: bool = False,

        # Optional alternative dataset
        load_marta_s: bool = False,
) -> xr.Dataset:
    """
    Loads the Dragonera insitu dataset using xarray.

    Parameters
    ----------
    time_selector: str | slice[str] | None, optional
        Time selection applied using xarray's `.sel()`. Can be a string (e.g., `"1993-01-21"`) or a slice
        (e.g., `slice("1993-01-21", "1993-01-25")`). If None, no time filtering is applied.
    
    nightime_data: bool, default=False
        If `True`, only the data between 9 p.m. and 6 a.m. (excluded) will be kept in the dataset.
        This is meant to match the CMEMS-SST time criterion.
    
    daily_avg: bool, default=False
        If `True`, daily averages the dataset.
        If nan value is present into a day, the resulting average of that day will be nan.
    
    drop_depth: bool, default=False
        If `False`, the depth coordinate will be set as 3m. If `True`, the depth is discarded.
    
    load_marta_s: bool, default=False
        If `True`, the dataset of Marta is loaded instead of the PdE dataset.
    
    Returns
    ----------
    ds_dragonera_insitu: xr.Dataset
        The loaded Dragonera insitu dataset.
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

    # Adding attributes
    ds_dragonera_insitu.time.attrs["reference"] = "GMT"
    ds_dragonera_insitu.T.attrs["unit"] = "°C"

    # Printing the good news
    details = []
    if nightime_data: details.append("nightime")
    if daily_avg: details.append("daily averaged")
    details = (" (" + ', '.join(details) + ")") if details != [] else ''
    
    # Printing the good news
    print(f"Loaded Dragonera insitu dataset{details}.")

    # Returning the final dataset
    return ds_dragonera_insitu

def save_dataset_to_nc(
        # Dataset
        ds: xr.Dataset,
        file_path: str,

        # Options
        progress_bar: bool = True,
        profilers: bool = False,
) -> xr.Dataset:
    """
    Basic function to save a dataset using xarray.

    Parameters
    ----------
    ds: xr.Dataset
        Dataset to save.

    file_path: str 
        Filepath where the dataset should be saved.
    
    progress_bar: bool, default=True
        If `True`, shows a progress bar when saving the dataset.
    
    profilers: bool, default=True
        If `True`, shows the profilers after the dataset has been saved.
    
    Returns
    ----------
    ds: xr.Dataset
        The computed dataset.
    """

    # Compute the dataset using a progress bar
    if progress_bar:
        with ProgressBar():
            ds = ds.compute()
            ds.to_netcdf(file_path)
    
    # Compute the dataset using a progress bar and other profilers
    elif progress_bar and profilers:
        with (ProgressBar(), ResourceProfiler() as rprof, Profiler() as prof):
            ds = ds.compute()
            ds.to_netcdf(file_path)
        
        # Save resource profiles as html file
        from datetime import datetime
        path = os.path.join(codigos_location, ".dask_profiles", datetime.now().isoformat() + ".html")
        visualize([prof, rprof], path, show=True, save=True)

    # Compute the dataset using no visualisation tools
    else:
        ds = ds.compute()
        ds.to_netcdf(file_path)

    # Returning the final dataset
    return ds

def load_nc_to_dataset(
        file_path: str,
        name_dict: dict = None,
        remove_bay_of_biscay: bool = False,
        region_selector: str | None = None,
        **kwargs
) -> xr.Dataset:
    """
    Loads a MHWs dataset using xarray.

    Parameters
    ----------
    file_path: str 
        Filepath from where the dataset should be loaded.

    remove_bay_of_biscay: bool, default=True
        If `True`, all the values of the dataset in the bay of biscay is discarded.
        Note that the coordinates must be named `lon` and `lat` for it to work.
    
    Returns
    ----------
    ds: xr.Dataset
        The loaded dataset.
    
    Other Parameters
    ----------------
    **kwargs
        All other arguments are passed to `xr.open_dataset`
    """

    # Loading the dataset
    ds = xr.open_dataset(
        file_path,
        decode_times=True,
        **kwargs
    )

    if name_dict:
        ds = ds.rename(name_dict)

    if region_selector == "balears":
        lon_selector = slice(-0.9, 5.1)
        lat_selector = slice(37.6, 41.1)
        ds = ds.sel(lat=lat_selector, lon=lon_selector)

    # Remove bay of Biscay
    if remove_bay_of_biscay:
        ds = ds.where((ds.lon > 0) | (ds.lat < 42))
    
    # Returning the final dataset
    return ds