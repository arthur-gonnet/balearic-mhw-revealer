########################################################################################################################
########################################################################################################################

###### USER NOTES:
# This script gather all the codes to load and save datasets using xarray. 

## Functions description
#
# - load_cmems_sst(...):
#       Loads the CMEMS-SST dataset using xarray.
#
# - save_cmems_sst_trend(...):
#       Saves a CMEMS-SST trend dataset using xarray.
#
# - load_cmems_sst_trend(...):
#       Loads a CMEMS-SST trend dataset using xarray.
#
# - save_cmems_sst_clim(...):
#       Saves a CMEMS-SST climatology dataset using xarray.
#
# - load_cmems_sst_clim(...):
#       Loads a CMEMS-SST climatology dataset using xarray.

########################################################################################################################
##################################### IMPORTS ##########################################################################
########################################################################################################################

# Basic imports
import os
import glob as glob

# Advanced imports
import xarray as xr

# Local imports
from data_managers.load_save_basics import load_nc_to_dataset, save_dataset_to_nc

# Load relative file paths
# These locations must be accurate to load the data
codigos_location = os.path.abspath(os.getcwd())
datos_location = os.path.join(codigos_location, "..", "Datos")

########################################################################################################################
##################################### USER INPUT #######################################################################
########################################################################################################################


# USER INPUT - CMEM-SST paths
cmems_sst_folder_path   = os.path.join(datos_location, "CMEMS-SST", "MEDITERRANEAN", "SST-L4-REP-HR", "DATA-NEW", "DAILY")
cmems_sst_ds_pattern    = os.path.join(cmems_sst_folder_path, "{year}", "SST_MED_SST_L4_REP_OBSERVATIONS_010_021_y{year}m{month}.nc")

# USER INPUT - CMEM-SST processed datasets paths
cmems_sst_trend_path    = os.path.join(datos_location, "CMEMS-SST", "processed", "cmems_sst_trend_1982_2023.nc")
cmems_sst_clim_pattern  = os.path.join(datos_location, "CMEMS-SST", "processed", "cmems_sst_clim_{region}_{year_start}_{year_end}.nc")


########################################################################################################################
##################################### CODES ############################################################################
########################################################################################################################


def load_cmems_sst(
        # Time preselection
        years: list[int] = range(1982, 2024),
        months: list[int] = range(1, 13),

        # Selectors
        time_selector: str | slice | None = None,
        lon_selector: float | slice | None = None,
        lat_selector: float | slice | None = None,
        region_selector: str | None = None,
        
        # Dataset options
        chunks: int | dict | str | None = 'auto',
        only_sst: bool = True,
) -> xr.Dataset:
    """
    Loads the CMEMS-SST dataset using xarray. By default, the entire dataset is loaded,
    but optional spatio-temporal subsettings are available for faster loading times.

    Parameters
    ----------
    years: list[int | str], default=['*']
        Years to load from the dataset, as integers (e.g., `range(1983, 1985)` or `[1983]`).
        Use `['*']` to load all available years (1982-2023) as it uses a glob pattern.
        This parameter defines the files to load, so loading time is highly depends on this range size.
    
    months: list[int | str], default=['*']
        Months to load from the dataset, as integers (e.g., `range(1, 5)` or `[1]`).
        Use `['*']` to load all months (1-12) as it uses a glob pattern.
        This parameter defines the files to load, so loading time is highly depends on this range size.
    
    time_selector: str | slice[str] | None, optional
        Time selection applied using xarray's `.sel()`. Can be a string (e.g., `"1993-01-21"`) or a slice
        (e.g., `slice("1993-01-21", "1993-01-25")`). If `None`, no time filtering is applied.
    
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
    
    only_sst: bool, default=True
        If `True`, all other variables than SST will be discarded.
    
    Returns
    ----------
    ds_cmems_sst: xr.Dataset
        The CMEMS-SST dataset.
    """

    # Selecting the files to load
    files: list[str] = []

    # Using temporal preslection of datasets
    for year in years:
        if year != '*' and (int(year) < 1982 or int(year) > 2023):
            print(f"Incorrect year selection for CMEMS-SST dataset, got {years}, expected a year range between 1982 and 2023")
        
        for month in months:
            if month != '*' and (int(month) < 1 or int(month) > 12):
                print(f"Incorrect month selection for CMEMS-SST dataset, got {months}, expected a month range between 1 and 12")
            
            month_str: str = str(month) if (month == '*' or month > 9) else ('0'+str(month))
            pattern: str = cmems_sst_ds_pattern.format(year=year, month=month_str)
            files.extend(glob.glob(pattern))
    
    # Load region selection
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
    Loads a CMEMS-SST climatology dataset using xarray.

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