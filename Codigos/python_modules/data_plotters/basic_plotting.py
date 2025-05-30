

import numpy as np
import xarray as xr

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from python_modules.utils import nice_range, override_value


def subplot(
        subplots_shape: tuple[int, int],
        subplots_settings: list[dict],

        figsize: tuple[float, float] = (10,6),
        fig_fontsize: int = 14,
        fig_title: str | None = None,

        fig_cbar: bool = False,
        fig_vmin: float | None = None,
        fig_vmax: float | None = None,
        fig_cmap: str | None = None,
        fig_cbar_unit: str | None = None,
        fig_cbar_pad=0.005,
        fig_cbar_fraction=0.025,

        # Parameter for saving the figure
        save_plot: bool = False,
        save_path: str = "",
        dpi: int = 160,

        show_plots: bool = False,
):
    if fig_cbar:
        if fig_vmin is None:
            fig_vmin = np.nanmin([np.nanmin(subplot_setting['variable_data']) for subplot_setting in subplots_settings])
        
        if fig_vmax is None:
            fig_vmax = np.nanmax([np.nanmax(subplot_setting['variable_data']) for subplot_setting in subplots_settings])
            fig_vmin, fig_vmax = nice_range(fig_vmin, fig_vmax)

    with plt.rc_context({'font.size': fig_fontsize}):
        fig = plt.figure(figsize=figsize, constrained_layout=True)

        # fig.subplots_adjust(wspace=0.4, hspace=0.4)

        for subplot_setting in subplots_settings:
            pos = subplot_setting.pop('pos')
            plotting_func = subplot_setting.pop('func')

            projection_mode = ccrs.PlateCarree()
            ax = fig.add_subplot(subplots_shape[0], subplots_shape[1], pos, projection=projection_mode)

            override_value(subplot_setting, "vmin", fig_vmin)
            override_value(subplot_setting, "vmax", fig_vmax)
            override_value(subplot_setting, "cmap", fig_cmap)
            override_value(subplot_setting, "fontsize", fig_fontsize)

            if fig_cbar:
                override_value(subplot_setting, "add_cbar", False)

            plotting_func(fig=fig, ax=ax, **(subplot_setting))

        if not fig_title is None:
            fig.suptitle(fig_title)

        if fig_cbar:
            norm = mcolors.Normalize(vmin=fig_vmin, vmax=fig_vmax)
            mappable = cm.ScalarMappable(norm=norm, cmap=fig_cmap)
            fig.colorbar(mappable, label=fig_cbar_unit, ax=fig.axes, fraction=fig_cbar_fraction, pad=fig_cbar_pad)


        # If saving the figure is asked
        if save_plot:
            fig.savefig(save_path, format="png", dpi=dpi)
        
        if show_plots:
            plt.show()
            plt.clf()
            plt.close("all")

def plot_map(
    # Parameter of data
    lon, lat,
    variable_data,

    # Parameter of figure
    fig=None,figsize=(10,6),
    ax: plt.Axes = None,

    # Parameter of text
    fontsize: int = 14,
    title: str="",

    # Parameter of colorbar
    add_cbar: bool=True,
    cbar_unit: str="",
    cbar_orientation: str="vertical",
    cbar_shrink: float=1,
    cbar_pad: float = 0.005,
    cbar_inversed: bool=False,
    vmin=None, vmax=None,

    bottom_labels = True,
    left_labels = True,

    # Parameters of graph
    extent="balears",
    aspect: float=1.3,
    cmap="managua",
    norm="linear",
    contours_levels=None,
    contours_labels_fontsize=8,
    xticks=None,
    yticks=None,

    # Parameter for saving the figure
    save_plot: bool=False,
    save_path: str="",
    dpi: int=160,

    show_plots: bool = False,
):
    """Creates a map plot of a variable over an area extent

    To use it directly, you must pass a :class:`pandas.DataFrame` as returned by a :class:`argopy.DataFetcher.index` or :class:`argopy.IndexFetcher.index` property::

        from argopy import IndexFetcher
        df = IndexFetcher(src='gdac').region([-80,-30,20,50,'2021-01','2021-08']).index
        bar_plot(df, by='profiler')

    Example
    --------
        from utils.plot_utils import plot_map

        ArgoSet = DataFetcher(mode='expert').float([6902771, 4903348]).load()
        ds = ArgoSet.data.argo.point2profile()
        df = ArgoSet.index

        scatter_map(df)
        scatter_map(ds)
        scatter_map(ds, hue='DATA_MODE')
        scatter_map(ds, hue='PSAL_QC')

    Parameters
    ----------
    df: :class:`pandas.DataFrame`
        As returned by a fetcher index property
    by: str, default='institution'
        The profile property to plot
    style: str, optional
        Define the Seaborn axes style: 'white', 'darkgrid', 'whitegrid', 'dark', 'ticks'

    Returns
    -------
    fig: :class:`matplotlib.figure.Figure`
    ax: :class:`matplotlib.axes.Axes`

    Other Parameters
    ----------------
    markersize: int, default=36
        Size of the marker used for profiles location.
    markeredgesize: float, default=0.5
        Size of the marker edge used for profiles location.
    markeredgecolor: str, default='default'
        Color to use for the markers edge. The default color is 'DARKBLUE' from :class:`argopy.plot.ArgoColors.COLORS`

    # kwargs
    #     All other arguments are passed to :class:`matplotlib.figure.Figure.subplots`
    """

    with plt.rc_context({'font.size': fontsize}):
        projection_mode = ccrs.PlateCarree()

        # If the figure is not configured, initialise a new one
        if fig is None or ax is None:
            fig = plt.figure(figsize=figsize, layout='tight')
            ax = fig.add_subplot(111, projection=projection_mode)
        
        # Changing number of ticks
        # xticks = ax.get_xticks()
        # ax.set_xticks(xticks[::len(xticks) // 4]) # set new tick positions
        # ax.tick_params(axis='x', rotation=30) # set tick rotation
        # ax.margins(x=0) # set tight margins
        # ax.locator_params('both', nbins=2)
        # ax.set_xticks(xticks)

        # Change the grid settings and the coordinates labels
        gl = ax.gridlines(linestyle=':', draw_labels=True, alpha=0.5)
        gl.top_labels = False
        gl.right_labels = False
        gl.bottom_labels = bottom_labels
        gl.left_labels = left_labels
        
        if xticks:
            xlocator = MaxNLocator(nbins=xticks)
            gl.xlocator = xlocator
        
        if yticks:
            ylocator = MaxNLocator(nbins=yticks)
            gl.ylocator = ylocator
        
        # gl.xlabel_style = {'size': 55, 'color': 'black'}
        # gl.ylabel_style = {'size': 55, 'color': 'black'}
        # ax.grid()


        if extent == "med":
            extent = [-6, 36.33, 30, 46]
        elif extent == "balears":
            extent = [-1, 5, 37.7, 41]

        # Choosing the plot extent and aspect
        if not extent is None:
            ax.set_extent(extent, crs=projection_mode)
        if not aspect is None:
            ax.set_aspect(aspect) # This is a trick for performance sake, giving faster results, it imitates Mercator projection in 10sec vs 60sec


        # Plot the variable
        pcm_opts = {
            "transform": projection_mode,
            "cmap": cmap,
            "norm": norm,
            
            # Extra
            # "shading": "nearest"
        }

        if vmin is None:
            vmin = None if np.isnan(variable_data).all() else np.nanmin(variable_data)
        
        if vmax is None:
            vmax = None if np.isnan(variable_data).all() else np.nanmax(variable_data)
            vmin, vmax = nice_range(vmin, vmax)

        pcm = ax.pcolormesh(lon, lat, variable_data, zorder=-1, shading="auto", vmin=vmin, vmax=vmax, **pcm_opts)
        ax.set_title(title)

        # Shows the isolines on the plots
        if not contours_levels is None:
            if isinstance(contours_levels, int):
                if vmin is None:
                    vmin = np.nanmin(variable_data)
                if vmax is None:
                    vmax = np.nanmax(variable_data)

                # If the contour level is a int, generate nice contour levels using MaxNLocator
                locator = MaxNLocator(contours_levels)
                contours_levels = locator.tick_values(vmin, vmax)
            
            contours = ax.contour(lon, lat, variable_data, levels=contours_levels, colors='k', linewidths=1, alpha=0.6, transform=projection_mode, zorder=-1)
            ax.clabel(contours, levels=contours_levels, colors='k', fontsize=contours_labels_fontsize, zorder=-1)


        # Add colorbar
        if add_cbar:
            cbar = fig.colorbar(pcm, ax=ax, label=cbar_unit, orientation=cbar_orientation, shrink=cbar_shrink, pad=cbar_pad)

            if cbar_inversed:
                cbar.ax.invert_yaxis()


        # Additions of coastlines and land
        ax.coastlines(resolution="10m")# resolution="10m", color="k")       # Coastlines
        ax.add_feature(cfeature.LAND)# color="gray")        # Land
        
        # If saving the figure is asked
        if save_plot:
            fig.savefig(save_path, format="png", dpi=dpi)
        
        if show_plots:
            plt.show()
            plt.clf()
            plt.close("all")

def plot_timeserie(
        # Datasets options
        vars: dict[str, xr.DataArray] | xr.DataArray,
        times: dict[str, xr.DataArray] | xr.DataArray,
        vars_stds: dict[str, xr.DataArray] | xr.DataArray | None = None,
        labels: dict[str, str] | str | None = 'auto',
        colors: dict[str, str] | str | None = None,
        nans_to_zero: bool = False,

        # Figure options
        fig: plt.Figure | None = None,
        ax: plt.Axes | None = None,
        figsize: tuple[int, int] = (18, 5),

        # Axe text options
        fontsize: int = 14,
        title: str | None = None,
        unit: str | None = None,

        # Axe options
        grid: bool = True,
        xlim: tuple[int|None, int|None] = (None, None),
        ylim: tuple[int|None, int|None] = (None, None),

        # Saving options
        show_plot: bool = False,
        save_plot: bool = False,
        save_path: str = "",
):
    """
    Plots a time serie with the defined settings.

    Parameters
    ----------
    vars: dict[str, xr.DataArray]
        Dictionary attributing to a name a data array representing the variable to plot.
    
    times: dict[str, xr.DataArray]
        Dictionary attributing to a name a data array representing the time array associated with the variable to plot.
    
    vars_stds: dict[str, xr.DataArray] | None, default=None
        Dictionary attributing to a name a data array representing the std array associated with the variable to plot.
    
    labels: dict[str, str] | None, default=None
        Dictionary attributing to a name the label to assign to the variable to plot.
    
    colors: dict[str, str] | None, default=None
        Dictionary attributing to a name the color to assign to the variable to plot.
    
        
    fig: plt.Figure | None, default=None
        Figure to use if previously created. None value will create a new figure only if ax is None.
        If ax is not None, no new figure will be created.
    
    ax: plt.Axes | None, default=None
        Axes to use if previously created. None value will create a new axes.
    
    figsize: tuple[int, int], default=(18, 5)
        Figure size to use if creating a new figure.
    
    
    fontsize: int = 14
        Font size to be used in the plot.

    title: str|None = None,
    unit: str|None = None,

    # Axe options
    grid: bool = True,
    xlim: tuple = (None, None),
    ylim: tuple = (None, None),

    # Saving options
    save_plot: bool = False,
    save_path: str = "",

    years: list[int|str], default=["*"], optional
        Years to load from the dataset, as integers or strings (e.g., range(1983, 1985) or ["1983"]).
        Use ["*"] to load all available years (1982-2023) as it uses a glob pattern.
    
    months: list[int|str], default=["*"], optional
        Months to load from the dataset, as integers or strings (e.g., range(1, 5) or ["1"]).
        Use ["*"] to load all months (1-12) as it uses a glob pattern.
    
    time_selector: str | slice[str] | None, default=None, optional
        Time selection applied using xarray's `.sel()`. Can be a string (e.g., "1993-01-21") or a slice
        (e.g., slice("1993-01-21", "1993-01-25")). If None, no time filtering is applied.
    
    lon_selector: float | slice[float] | None, default=None, optional
        Longitude selector applied using xarray's `.sel()`. Accepts a float or a slice.
    
    lat_selector: float | slice[float] | None, default=None, optional
        Latitude selector applied using xarray's `.sel()`. Accepts a float or a slice.
    
    only_sst: bool, default=False, optional
        If True, all other variables than SST will be discarded.
    
    Returns
    ----------
    ds_cmems_sst: xr.Dataset
        The loaded CMEMS-SST dataset with optional spatio-temporal subsetting.
    """

    # min_date = np.min([ds_vars["time"][0] for ds_vars in datasets_vars])
    # max_date = np.min([ds_vars["time"][-1] for ds_vars in datasets_vars])

    with plt.rc_context({'font.size': fontsize}):
        # Plotting
        plt.figure(figsize=figsize, layout='tight')
        
        # Should handle single dataset plotting
        if isinstance(vars, dict):
            first_var = vars[next(iter(vars))]

            if unit == None and isinstance(first_var, xr.DataArray) and first_var.attrs.get("unit"):
                unit = first_var.attrs.get("unit")
            
            for var in vars:
                opts_args = {}
                if colors: opts_args["color"] = colors[var]
                
                # Labelling
                if labels == 'auto':
                    opts_args["label"] = var
                elif labels:
                    opts_args["label"] = labels[var]
                
                if nans_to_zero:
                    vars[var] = np.nan_to_num(vars[var])
                
                plt.plot(times[var], vars[var], lw=1, **opts_args)

                if xlim == (None, None):
                    xlim = (
                        np.nanmin([np.nanmin(times[name]) for name in times]),
                        np.nanmax([np.nanmax(times[name]) for name in times])
                    )
            
        else:
            if unit == None and isinstance(vars, xr.DataArray) and vars.attrs.get("unit"):
                unit = vars.attrs.get("unit")

            if labels == 'auto': labels=None

            if nans_to_zero:
                vars = np.nan_to_num(vars)

            plt.plot(times, vars, color=colors, label=labels, lw=1)  

            if xlim == (None, None):
                xlim = (
                    np.nanmin(times),
                    np.nanmax(times)
                )      

        plt.title(title)
        plt.ylabel(f"[{unit}]")

        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.grid(alpha=0.2)
        if labels: plt.legend()

        # plt.tight_layout()

        if show_plot:
            plt.show()
            plt.clf()
            plt.close("all")