########################################################################################################################
########################################################################################################################

###### USER NOTES:
# This script gather all the codes to load and save datasets using xarray. 

## Functions description
#
# - save_mhws_dataset(...):
#       Saves a MHWs dataset using xarray.
#
# - load_mhws_dataset(...):
#       Loads a MHWs dataset using xarray.

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


# USER INPUT - MHWs datasets pathes
mhws_dataset_pattern = os.path.join(datos_location, "mhws", "{dataset}_mhws_{region}_{clim_start}_{clim_end}.nc")


########################################################################################################################
##################################### CODES ############################################################################
########################################################################################################################


def save_mhws_dataset(
        # Dataset
        ds_mhws: xr.Dataset,

        # File path parameters
        dataset_used: str,
        region: str = "balears",
        clim_period: tuple[int, int] = (1987, 2021),

        # Options
        progress_bar: bool = True,
        profilers: bool = False
) -> xr.Dataset:
    """
    Saves a MHWs dataset using xarray.

    Parameters
    ----------
    ds_mhws: xr.Dataset
        MHWs dataset to save.

    dataset_used: str 
        Describes the dataset from which the MHWs computations were performed.
        Can be `"cmems_sst"`, `"med_mfc_bot"` or `"med_mfc_50m"` for example.

    region: str, default="balears"
        Region in which the MHWs computations were performed. Can be `"balears"` or `"med"` for example.

    clim_period: tuple[int, int], default=(1987,2021)
        Climatology period used for the MHWs computations.
    
    progress_bar: bool, default=True
        If True, shows a progress bar when saving the dataset.
    
    profilers: bool, default=True
        If True, shows the profilers after the dataset has been saved.
    
    Returns
    ----------
    ds_mhws: xr.Dataset
        The computed MHWs dataset.
    """
    
    # Getting dataset's associated filepath
    mhws_dataset_path = mhws_dataset_pattern.format(
        dataset = dataset_used,
        region = region,
        clim_start = clim_period[0],
        clim_end = clim_period[1]
    )

    # Printing the good news
    print(f"Saving MHWs dataset to {mhws_dataset_path}")

    # Saving dataset to nc
    ds_mhws = save_dataset_to_nc(
        ds           = ds_mhws,
        file_path    = mhws_dataset_path,

        progress_bar = progress_bar,
        profilers    = profilers
    )

    # Printing the good news
    print(f" -> Saved!")

    # Returning the dataset after Dask calculations have been performed
    return ds_mhws

def load_mhws_dataset(
        # File path parameters
        dataset_used: str,
        region: str = "balears",
        clim_period: tuple[int, int] = (1987, 2021),
) -> xr.Dataset:
    """
    Loads a MHWs dataset using xarray.

    Parameters
    ----------
    dataset_used: str 
        Describes the dataset from which the MHWs computations were performed.
        Can be `"cmems_sst"`, `"med_mfc_bot"` or `"med_mfc_50m"` for example.

    region: str, default="balears"
        Region in which the MHWs computations were performed. Can be `"balears"` or `"med"` for example.

    clim_period: tuple[int, int], default=(1987,2021)
        Climatology period used for the MHWs computations.
    
    Returns
    ----------
    ds_mhws: xr.Dataset
        The loaded MHWs dataset.
    """

    # Getting dataset's associated filepath
    mhws_dataset_path = mhws_dataset_pattern.format(
        dataset = dataset_used,
        region = region,
        clim_start = clim_period[0],
        clim_end = clim_period[1]
    )

    # Loading the dataset
    ds_mhws = load_nc_to_dataset(
        file_path = mhws_dataset_path,
        remove_bay_of_biscay = True
    )

    # Printing the good news
    print("Loaded MHWs dataset.")
    
    # Returning the final dataset
    return ds_mhws
