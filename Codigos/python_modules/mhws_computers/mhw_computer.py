########################################################################################################################
########################################################################################################################

###### USER NOTES:
# This script gathers all the codes to load and save datasets using xarray. 

## Functions description
#
# - compute_mhw_maps_apply_ufunc(...):
#       ...
#
# - compute_mhw_wrapped(...):
#       ...

########################################################################################################################
##################################### IMPORTS ##########################################################################
########################################################################################################################

# Basic imports
import os
import glob as glob

# Advanced imports
import numpy as np
import xarray as xr

# Local imports
import mhws_computers.marineHeatWaves as mhw
import options as opts

########################################################################################################################
##################################### CODES ############################################################################
########################################################################################################################


def compute_mhw_yearly(
        ds: xr.Dataset | xr.DataArray,
        using_dataset: str,
        var_name: str = "T",

        clim_period: tuple[int, int] = (1987,2021),
):
    if isinstance(ds, xr.Dataset):
        ds = ds[var_name]

    stacked = False

    # Stacking lon and lat for perfomance sake
    if ds.lon.size > 1 and ds.lat.size > 1:
        ds = ds.stack(pos=("lon", "lat"))
        stacked = True
    
    # Years of the dataset
    years = np.unique(ds["time.year"].values)

    # Computing the mhws
    #   result contains the resulting mhws statistics
    #   shape is (mhws_stats=18, pos, years)
    results = xr.apply_ufunc(
        # Inputs
        compute_mhw_yearly_wrapped,
        ds.time,
        ds.chunk({dim: -1 if dim == "time" else "auto" for dim in ds.dims}),
        kwargs={
            "clim_period": clim_period
        },

        # Dimensions of input and output
        input_core_dims     = [['time'], ['time']],
        output_core_dims    = [(['year'])] * len(opts.mhws_stats),
        
        # Type and size of output
        output_dtypes       = [float] * len(opts.mhws_stats),
        dask_gufunc_kwargs  = dict(
            output_sizes = {
                'year': len(years)
            }
        ),

        # Dask Options
        vectorize   = True,
        dask        = 'parallelized',
    )

    # print(results)

    # Assigning coordinates
    output_ds = xr.Dataset({
        var: data.assign_coords(year=years)
        for var, data in zip(opts.mhws_stats, results)
    })

    # Unstacking pos
    if stacked:
        output_ds = output_ds.unstack("pos")

    # Changing dimensions order
    output_ds = output_ds.transpose('lat', 'lon', 'year')
    # output_ds = output_ds.transpose('year')

    for stat in opts.mhws_stats:
        output_ds[stat].attrs["shortname"]  = opts.mhws_stats_shortname[stat]
        output_ds[stat].attrs["longname"]   = opts.mhws_stats_longname[stat]
        output_ds[stat].attrs["unit"]       = opts.mhws_stats_units[stat]

    # Adding attributes
    output_ds.attrs['climatologyPeriod'] = f'{clim_period[0]}-{clim_period[1]}'
    output_ds.attrs['description'] = opts.mhw_dataset_description

    match using_dataset.lower():
        case "cmems_sst":
            output_ds.attrs['acknowledgment'] = opts.cmems_sst_acknowledgment

        case "cmems_mfc":
            output_ds.attrs['acknowledgment'] = opts.cmems_mfc_acknowledgment

    return output_ds

def compute_mhw_yearly_wrapped(t: np.array, sst: np.array, clim_period: tuple[int, int]):
    # print(f"received time array of type {type(t)}, shape {t.shape}")
    # print(f"received sst array of type {type(sst)}, shape {sst.shape}")

    # t = ((t - np.datetime64("0001-01-01")) / np.timedelta64(1, "D")).astype(int)
    t = t.astype('datetime64[D]').astype(int) + 719163
    sst = sst.copy()

    mhws, clim = mhw.detect(t, sst, climatologyPeriod=clim_period, splitYearOverlappingMHWs=True)
    mhwBlock = mhw.blockAverage(t, mhws, clim, temp=sst)

    return tuple(mhwBlock[stat] for stat in opts.mhws_stats)


    ## All events approach (not yearly stats)

mhws_all_events_stats = [
    'time_start', 'time_end', 'time_peak',
    # 'date_start', 'date_end', 'date_peak',
    'index_start', 'index_end',
    'index_peak', 'duration', 'duration_moderate', 'duration_strong', 'duration_severe', 'duration_extreme',
    'intensity_max', 'intensity_mean', 'intensity_var', 'intensity_cumulative', 'intensity_max_relThresh',
    'intensity_mean_relThresh', 'intensity_var_relThresh', 'intensity_cumulative_relThresh', 'intensity_max_abs',
    'intensity_mean_abs', 'intensity_var_abs', 'intensity_cumulative_abs', 'category', 'rate_onset', 'rate_decline',]
clim_keys = ['thresh', 'seas', 'missing']
    
def compute_mhw_all_events(
        ds: xr.Dataset | xr.DataArray,
        using_dataset: str,
        var_name: str = "T",

        clim_period: tuple[int, int] = (1987,2021),
):
    if isinstance(ds, xr.Dataset):
        ds = ds[var_name]

    # dims_to_stack = [dim for dim in ds.dims if dim != "time"]
    # should_stack = sum([ds[dim].size > 1 for dim in dims_to_stack]) > 1
    stack = ds.lon.size > 1 and ds.lat.size > 1

    # Stacking lon and lat for perfomance sake
    if stack:
        ds = ds.stack(pos=("lon", "lat"))
    
    from dask.diagnostics import ProgressBar
    # As operations are not dask-oriented, need to compute
    with ProgressBar():
        ds = ds.compute()

    # Computing the mhws
    #   result contains the resulting mhws statistics
    #   shape is (things=33, pos, event_number/time)
    results = xr.apply_ufunc(
        # Inputs
        compute_mhw_all_events_wrapped,
        ds.time,
        ds,
        kwargs={
            'clim_period': clim_period
        },

        # Dimensions of input and output
        input_core_dims     = [['time'], ['time']],
        output_core_dims    = [['event_number']] * len(mhws_all_events_stats) + [['time']] * len(clim_keys),
        
        # Type and size of output
        output_dtypes       = [float] * (len(mhws_all_events_stats) + len(clim_keys)),
        # dask_gufunc_kwargs  = dict(
        #     output_sizes = {
        #         'event_number': ...,
        #         'time': ds.time.size,
        #     }
        # ),

        # Dask Options
        vectorize   = True,
        # dask        = 'parallelized', # Cannot as we don't know the events number
    )

    # print(results)

    # Assigning coordinates
    output_ds = xr.Dataset({
        var: data.assign_coords(event_number=range(len(results[0])))
        for var, data in zip(mhws_all_events_stats, results[:len(mhws_all_events_stats)])
    } | {
        'clim_'+var: data.assign_coords(time=ds.time)
        for var, data in zip(clim_keys, results[len(mhws_all_events_stats):])
    })

    # Unstacking pos
    if stack:
        output_ds = output_ds.unstack("pos")

    # Changing dimensions order
    # output_ds = output_ds.transpose('lat', 'lon', 'year')
    # output_ds = output_ds.transpose('year')

    # for stat in opts.mhws_stats:
    #     output_ds[stat].attrs["shortname"]  = opts.mhws_stats_shortname[stat]
    #     output_ds[stat].attrs["longname"]   = opts.mhws_stats_longname[stat]
    #     output_ds[stat].attrs["unit"]       = opts.mhws_stats_units[stat]

    # Adding attributes
    output_ds.attrs['climatologyPeriod'] = f'{clim_period[0]}-{clim_period[1]}'
    output_ds.attrs['description'] = opts.mhw_dataset_description

    match using_dataset.lower():
        case "cmems_sst":
            output_ds.attrs['acknowledgment'] = opts.cmems_sst_acknowledgment

        case "cmems_mfc":
            output_ds.attrs['acknowledgment'] = opts.cmems_mfc_acknowledgment

    return output_ds

def compute_mhw_all_events_wrapped(t: np.array, sst: np.array, clim_period: tuple[int, int]):
    # print(f"received time array of type {type(t)}, shape {t.shape}")
    # print(f"received sst array of type {type(sst)}, shape {sst.shape}")

    # t = ((t - np.datetime64("0001-01-01")) / np.timedelta64(1, "D")).astype(int)
    t = t.astype('datetime64[D]').astype(int) + 719163
    sst = sst.copy()

    print(t.shape, t)
    print(sst.shape, sst)
    mhws, clim = mhw.detect(t, sst, climatologyPeriod=clim_period)

    for stat in mhws_all_events_stats:
        print(f"{stat}: {'list' if isinstance(np.array(mhws[stat]), list) else np.array(mhws[stat]).dtype}, {np.array(mhws[stat]).shape}" )
    
    for key in clim_keys:
        print(f"{key}: {clim[key].dtype}, {clim[key].shape}" )
        
    categories = {'Moderate': 1, 'Strong': 2, 'Severe': 3, 'Extreme': 4}
    mhws["category"] = [categories[cat] for cat in mhws["category"]]

    print([mhws[stat] for stat in mhws_all_events_stats] + [clim[key] for key in clim_keys])
    return tuple([np.array(mhws[stat]) for stat in mhws_all_events_stats] + [clim[key].astype(float) for key in clim_keys])


def compute_mhw_ts(
        ds: xr.Dataset | xr.DataArray,
        clim_period: tuple[int, int] = (1987,2021),
        var_name: str = "T",

        blockAveraging: bool = False,
) -> tuple[dict, dict] | dict:

    if isinstance(ds, xr.Dataset):
        ds = ds[var_name]

    t = ds.time.values.astype('datetime64[D]').astype(int) + 719163
    sst = ds.values

    mhws, clim = mhw.detect(t, sst, climatologyPeriod=clim_period)

    if blockAveraging:
        mhwBlock = mhw.blockAverage(t, mhws, clim, temp=sst)
        return mhwBlock

    return mhws, clim





































""" (https://discourse.pangeo.io/t/question-on-dask-efficiency/2366/21?page=2)
# Create fake input data
lat_size, long_size = 100, 100
data = da.random.random_integers(0, 30, size=(1_000, long_size, lat_size), chunks=(-1, 10, 10))  # size = (time, longitude, latitude)
time = np.arange(730_000, 731_000)  # time in ordinal days

# define a wrapper for the climatology
def func1d_climatology(arr, time):
   _, point_clim = mhw.detect(time, arr)
   # return climatology
   return point_clim['seas']

# define a wrapper for the threshold
def func1d_threshold(arr, time):
   _, point_clim = mhw.detect(time, arr)
   # return threshold
   return point_clim['thresh']

# output arrays
full_climatology = da.zeros_like(data)
full_threshold = da.zeros_like(data)

climatology = da.apply_along_axis(func1d_climatology, 0, data, time=time, dtype=data.dtype, shape=(1000,))
threshold = da.apply_along_axis(func1d_threshold, 0, data, time=time, dtype=data.dtype, shape=(1000,))
climatology = climatology.compute()
threshold = threshold.compute()
"""

def compute_mhw_ds(ds: xr.Dataset | xr.DataArray, var_name: str="T", **kwargs):
    if isinstance(ds, xr.Dataset):
        ds = ds[var_name]
    
    # ds = ds.assign_coords({
    #     # "ordinal_time": ("time", ((ds.time - np.datetime64("0001-01-01")) / np.timedelta64(1, "D")).astype(int))
    #     # "ordinal_time": ("time", ((ds.time - np.datetime64("0001-01-01")) / np.timedelta64(1, "D")).astype(int))
    #     "ordinal_time": ("time", ds.time.values.astype('datetime64[D]').astype(int) + 719163)
    # })
    
    # return compute_mhw(ds.ordinal_time.values, ds.values, **kwargs)
    # return compute_mhw(ds.time.values.astype('datetime64[D]').astype(int) + 719163, ds.values, **kwargs)
    return compute_mhw(((ds.time.astype('datetime64[D]') - np.datetime64("0001-01-01")) / np.timedelta64(1, "D")).astype(int).values, ds.values, **kwargs)


def compute_mhw(
        t, sst,
        lon=np.nan, lat=np.nan, depth=np.nan,
        ds_attrs: dict={
            "description": "MHWs yearly statistics computed using the marineHeatWaves module for python developped by Eric C. J. Oliver."
        },
        blocks: bool = True,
        **kwargs
):
    mhws, clim = mhw.detect(t, sst, **kwargs)

    if not blocks:
        return mhws, clim
    
    mhwBlock = mhw.blockAverage(t, mhws, clim)

    datavars: dict = {}

    for stat in opts.stat_verbose:
        datavars[stat] = (
            ("lon", "lat", "depth", "year"),
            np.expand_dims(
                np.expand_dims(
                    np.expand_dims(
                        mhwBlock[stat],
                    axis=0),
                axis=0),
            axis=0),
            {
                "longname": opts.stat_verbose[stat],
                "unit": opts.stat_units[stat],
            }
        )
    
    ds = xr.Dataset(
        data_vars = datavars,
        coords = dict(
            lon=[lon],
            lat=[lat],
            depth=[depth],
            year=mhwBlock["years_start"].astype(int),
        ),
        attrs = ds_attrs,
    )

    return ds


def nan_mhwBlock_ds(lon, lat, depth, years, **kwargs):
    datavars: dict = {}

    for stat in opts.stat_verbose:
        datavars[stat] = (
            ("lon", "lat", "depth", "year"),
            np.expand_dims(
                np.expand_dims(
                    np.expand_dims(
                        np.full(len(years), np.nan),
                        # np.empty(len(years)),
                    axis=0),
                axis=0),
            axis=0),
            {
                "name": opts.stat_verbose[stat],
                "unit": opts.stat_units[stat],
            }
        )
    
    ds = xr.Dataset(
        data_vars = datavars,
        coords = dict(
            lon=[lon],
            lat=[lat],
            depth=[depth],
            year=years,
        ),
        attrs = dict(description="MHWs yearly statistics computed using the marineHeatWaves module for python developped by Eric C. J. Oliver."),
    )

    return ds
