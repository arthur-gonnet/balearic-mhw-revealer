




import math

import xarray as xr

import cmocean.cm as cmo

from data_plotters.basic_plotting import subplot, plot_map
from utils import not_null
import options as opts




def plot_1year_stats_mhws_maps(
        ds_mhws: xr.Dataset,
        year: int,
        stat_set: list[str] | str | int,
        ncols: int | None = None,
        custom_vlims: dict = {},
        custom_cmaps: dict = {},
        title = "MHWs statistics in {year}"
):
    """"
    stat_set: can be
            "basic", "advanced", "categories", "temp_stats", "all"
            ["...", "..."]
    """

    match stat_set:
        case "basic" | 1:
            stat_set = opts.mhws_basic_stats
            ncols = not_null(ncols, 2)
        case "advanced" | 2:
            stat_set = opts.mhws_advanced_stats
            ncols = not_null(ncols, 2)
        case "categories" | 3:
            stat_set = opts.mhws_categories_stats
            ncols = not_null(ncols, 2)
        case "temp_stats" | 4:
            stat_set = opts.mhws_temp_stats
            ncols = not_null(ncols, 2)
        case "all" | 5:
            stat_set = opts.mhws_stats
            ncols = not_null(ncols, 3)
        case _:
            ncols = not_null(ncols, 2)
    
    if isinstance(stat_set, str):
        stat_set = [stat_set]
    
    nrows = math.ceil(len(stat_set) / ncols)

    subplots_settings = []

    for i, stat in enumerate(stat_set):
        row = i // ncols
        col = i % ncols

        subplots_settings.append(dict(
            pos = i+1,
            func = plot_map,
            lon = ds_mhws.sel(year=year).lon,
            lat = ds_mhws.sel(year=year).lat,
            variable_data = ds_mhws.sel(year=year)[stat],

            # Parameter of text
            # fontsize = 22,
            title = opts.mhws_stats_shortname[stat],

            # Parameter of colorbar
            add_cbar = True,
            cbar_pad = 0.01,
            cbar_unit = f"[{opts.mhws_stats_units[stat]}]",
            vmin = custom_vlims.get(stat, [None])[0],
            vmax = custom_vlims.get(stat, [None,None])[1],
            cmap = custom_cmaps.get(stat, opts.mhws_stats_cmaps[stat]),

            left_labels = (col == 0),
            bottom_labels = (row == nrows - 1),
            xticks = 6,
            yticks = 4,
        ))

    subplot(
        (nrows, ncols),
        subplots_settings,
        
        figsize = (ncols*7.3, nrows*4.67),
        fig_fontsize = 20,

        fig_title=title.format(year=year),
        show_plots = True,
    )

def plot_1stat_years_mhws_maps(
        ds_mhws,
        years,
        stat,
        ncols=2,
        custom_vlims={}
):
    nrows = math.ceil(len(years) / ncols)

    subplots_settings = []

    for i, year in enumerate(years):
        row = i // ncols
        col = i % ncols

        subplots_settings.append(dict(
            pos = i+1,
            func = plot_map,
            lon = ds_mhws.sel(year=year).lon,
            lat = ds_mhws.sel(year=year).lat,
            variable_data = ds_mhws.sel(year=year)[stat],

            # Parameter of text
            # fontsize = 22,
            title = year,

            left_labels = (col == 0),
            bottom_labels = (row == nrows - 1),
            xticks = 6,
            yticks = 4,
        ))

    subplot(
        (nrows, ncols),
        subplots_settings,
        
        figsize = (ncols*7.3, nrows*4.67),
        fig_fontsize = 20,

        fig_cbar=True,
        fig_cmap=opts.mhws_stats_cmaps[stat],
        fig_cbar_unit=opts.mhws_stats_units[stat],
        fig_cbar_fraction=0.01,
        fig_vmin = custom_vlims.get(stat, [None])[0],
        fig_vmax = custom_vlims.get(stat, [None,None])[1],

        fig_title=f"{opts.mhws_stats_shortname[stat]}",
        show_plots = True,
    )


def plot_mhws_timeseries(
        ds_mhw,
        stat_set,
        ncols=2,
):
    ...
