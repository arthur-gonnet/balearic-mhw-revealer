
# Basic libraries
import os, sys

# This makes sure that the current directory is set to be the root Codigos/ directory
# TO BE UPDATED FOR EACH SCRIPT DEPENDING OF FILE TREE
codigos_path = os.path.abspath(os.getcwd())
sys.path.append(codigos_path)

# So the Datos and Figuras locations are just one folder up from the Codigos path
datos_location = os.path.join(codigos_path, '..', 'Datos')
figuras_location = os.path.join(codigos_path, '..', 'Figuras')

# Advanced libraries
import numpy as np
import pandas as pd
import xarray as xr


# USER INPUT - The initial and final path to the files
in_dragonera_insitu_csv_path = os.path.join(datos_location,"puerto_del_estado","raw_files","20623_38448_2820_WATER_TEMP_20061129084612_20250430074612.csv")
out_dragonera_insitu_nc_path = os.path.join(datos_location,"puerto_del_estado","dragonera_hist_temp.nc")

in_dragonera_insitu_marta_s_txt_path = os.path.join(datos_location,"puerto_del_estado","raw_files","IR_TS_MO_6100430_TEMPb.txt")
out_dragonera_insitu_marta_s_nc_path = os.path.join(datos_location,"puerto_del_estado","dragonera_hist_temp_marta_s.nc")


def puerto_csv_to_nc(
        fill_with_nans: bool = True,
        debug: bool = False,
        to_nc: bool = False,
):
    # Reading the initial csv file
    df: pd.DataFrame = pd.read_csv(in_dragonera_insitu_csv_path, delimiter="\t", header=1)
    if debug:print("Initial data: ", df, "", sep="\n")

    # Changing the columns names
    df = df.rename(columns={
        "Fecha (GMT)": "time",
        "Temperatura del Agua(ºC)": "T"
    }, errors="raise")
    if debug:print("Changed columns names: ", df, "", sep="\n")

    # Setting the time from string to datetime format
    df["time"] = pd.to_datetime(df["time"], format="%Y %m %d %H")
    if debug:print("Changed time format: ", df, "", sep="\n")

    # Setting the columns to represent coordinates
    df = df.set_index("time")
    if debug:print("Added time index: ", df, "", sep="\n")

    # Adding nan values to have a consistent time serie
    if fill_with_nans:
        full_index = pd.date_range(start=df.index.min(), end=df.index.max(), freq="1h")
        df = df.reindex(full_index)
        df.index.name = "time"

        # for date in pd.date_range(start=df.time.min(), end=df.time.max(), freq="1h"):
        #     if not date in df.time.unique():
        #         df = pd.concat([pd.DataFrame([[date, np.nan, np.nan]], columns=df.columns), df], ignore_index=True)
        if debug:print("Reindexed time with nans: ", df, "", sep="\n")

    # Converting to xarray dataset
    xrdf: xr.Dataset = df.to_xarray()
    if debug:print("Final xr dataset: ", xrdf, "", sep="\n")

    # Finally saving as netcdf file
    if to_nc:
        xrdf.to_netcdf(out_dragonera_insitu_nc_path)

    return df


def puerto_marta_s_txt_to_nc(
        fill_with_nans: bool = True,
        debug: bool = False,
        to_nc: bool = False,
):
    # Reading the initial csv file
    df: pd.DataFrame = pd.read_csv(
        in_dragonera_insitu_marta_s_txt_path,
        delimiter="  ",
        names=["Canal de obtencion de los datos", "time", "T", "one"],
        usecols=["time", "T"],
        engine="python"
    )
    if debug:print("Initial data: ", df, "", sep="\n")

    # Setting the time from string to datetime format
    df["time"] = pd.to_datetime(df["time"], format="%Y %m %d %H 00 00")
    if debug:print("Changed time format: ", df, "", sep="\n")

    # Setting the columns to represent coordinates
    df = df.set_index("time")
    if debug:print("Added time index: ", df, "", sep="\n")

    # Adding nan values to have a consistent time serie
    if fill_with_nans:
        full_index = pd.date_range(start=df.index.min(), end=df.index.max(), freq="1h")
        df = df.reindex(full_index)
        df.index.name = "time"

        # for date in pd.date_range(start=df.time.min(), end=df.time.max(), freq="1h"):
        #     if not date in df.time.unique():
        #         df = pd.concat([pd.DataFrame([[date, np.nan, np.nan]], columns=df.columns), df], ignore_index=True)
        if debug:print("Reindexed time with nans: ", df, "", sep="\n")

    # Converting to xarray dataset
    xrdf: xr.Dataset = df.to_xarray()
    if debug:print("Final xr dataset: ", xrdf, "", sep="\n")

    # Finally saving as netcdf file
    if to_nc:
        xrdf.to_netcdf(out_dragonera_insitu_marta_s_nc_path)

    return df

# def convert_csv_to_nc(
#         input_filepath: str,

#         kwargs_pd_read_csv: dict = {},
#         columns_mapper: dict = None,
#         handle_time: bool = False,

#         fill_with_nans: bool = True,
#         debug: bool = False,
#         to_nc: bool = False,
# ):
#     # Reading the initial csv file
#     df: pd.DataFrame = pd.read_csv(input_filepath, **kwargs_pd_read_csv)
#     if debug:print("Initial data: ", df, "", sep="\n")

#     # Changing the columns names
#     if not columns_mapper is None:
#         df = df.rename(columns=columns_mapper, errors="raise")
#         if debug:print("Changed columns names: ", df, "", sep="\n")

#     if handle_time:
#         # Setting the time from string to datetime format
#         df["time"] = pd.to_datetime(df["time"], format="%Y %m %d %H")
#         if debug:print("Changed time format: ", df, "", sep="\n")

#         # Setting the columns to represent coordinates
#         df = df.set_index("time")
#         if debug:print("Added time index: ", df, "", sep="\n")

#         # Adding nan values to have a consistent time serie
#         if fill_with_nans:
#             full_index = pd.date_range(start=df.index.min(), end=df.index.max(), freq="1h")
#             df = df.reindex(full_index)
#             df.index.name = "time"

#             # for date in pd.date_range(start=df.time.min(), end=df.time.max(), freq="1h"):
#             #     if not date in df.time.unique():
#             #         df = pd.concat([pd.DataFrame([[date, np.nan, np.nan]], columns=df.columns), df], ignore_index=True)
#             if debug:print("Reindexed time with nans: ", df, "", sep="\n")

#     # Converting to xarray dataset
#     xrdf: xr.Dataset = df.to_xarray()
#     if debug:print("Final xr dataset: ", xrdf, "", sep="\n")

#     # Finally saving as netcdf file
#     if to_nc:
#         xrdf.to_netcdf(out_dragonera_insitu_nc_path)

#     return df
