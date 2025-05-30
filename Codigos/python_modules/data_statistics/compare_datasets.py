

import numpy as np
import xarray as xr

from dask.diagnostics import ProgressBar

def basic_stats(
        da_1: xr.DataArray | xr.Dataset,
        da_2: xr.DataArray | xr.Dataset,
        var_name: str | None = None,
        precised_stats: bool = True,
):
    if isinstance(da_1, xr.Dataset):
        if var_name is None:
            print("ERROR. Differentiating dataset without specifying which variable to compare.")
            return
        else:
            da_1 = da_1[var_name]
            da_2 = da_2[var_name]
    
    with ProgressBar():
        da: xr.DataArray = (da_1 - da_2)
        npda: np.array = da.values

    if not da_1.attrs.get("unit"):
        print("No unit found for the DataArray, using °C.")

    unit = da_1.attrs.get('unit', '°C')

    stats = [
        {
            "name": "Bias",
            "value": np.nanmean(npda),
        },
        {
            "name": "RMSE",
            "value": np.nanmean(npda**2) ** 0.5,
        },
        {
            "name": "STDE",
            "value": np.nanstd(npda),
        }
    ]

    for stat in stats:
        if precised_stats:
            print(f"{stat['name']}: {stat['value']:.2f} {unit} ({stat['value']:.6f})")
        else:
            print(f"{stat['name']}: {stat['value']:.2f} {unit}")