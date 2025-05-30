
import os

import matplotlib.pyplot as plt

import numpy as np
import xarray as xr

from dask.diagnostics import ProgressBar, ResourceProfiler, Profiler, visualize

# from dask.distributed import Client
# client = Client()
codigos_location = os.path.abspath(os.getcwd())
datos_location = os.path.join(codigos_location, "..", "Datos")


def calculate_linear_trend(ds, var_name):

    # Converting time as float year unit
    ds["time"] = ds['time.year'].data + (ds['time.dayofyear'].data - 1) / 365.25

    # Operations to calculate the trend
    ds_trend = ds[var_name].polyfit(dim='time', deg=1)
    ds_trend = ds_trend.sel(degree=1, drop=True)
    ds_trend = ds_trend.rename({"polyfit_coefficients" : "sst_trend"})

    return ds_trend


def calculate_climatology(ds, time_name="time"):
    # Operations to calculate the climatology
    ds_clim = ds.groupby(f"{time_name}.dayofyear").mean()

    # Calculate 10th and 90th percentile
    

    return ds_clim