


from matplotlib.ticker import MaxNLocator


def nice_range(
    vmin: float | None,
    vmax: float | None,
    nbins: int = 10
) -> tuple[float, float]:
    """
    Uses matplotlib's MaxNLocator to compute nice rounded range and ticks.
    
    Parameters
    ----------
    vmin: float | None
        Initial vmin value. Note that if one of vmin or vmax is None, the range is returned unchanged.
    
    vmax: float | None
        Initial vmin value. Note that if one of vmin or vmax is None, the range is returned unchanged.
    
    nbins: int | 'auto', default: 10
        Matplotlib's MaxNLocator parameter.
        Maximum number of intervals; one less than max number of ticks.
        If the string 'auto', the number of bins will be automatically determined based on the length of the axis.

    Returns:
        nice_vmin, nice_vmax
    """

    if vmin is None or vmax is None:
        return vmin, vmax

    if vmin == 0 and vmax == 0:
    # if (vmin < 1e-10 and vmin > 1e-10 and vmax < 1e-10 and vmax > 1e-10) or (vmin is None and vmax is None):
        return 0., 1.
    
    
    locator = MaxNLocator(nbins=nbins, prune=None)
    ticks = locator.tick_values(vmin, vmax)

    return ticks[0], ticks[-1] # ticks[1] - ticks[0] is the step


def override_value(dict, key, value):
    if not value is None:
        dict[key] = value

def not_null(value, default):
    if value is None:
        return default
    return value