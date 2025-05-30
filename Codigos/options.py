########################################################################################################################
########################################################################################################################

###### USER NOTES:
# This script gathers some options used along the codes of this repository. 

########################################################################################################################
##################################### USER INPUT ##########################################################################
########################################################################################################################


    # Maps options

# CURRENTLY NOT USED IN CODE
extents = {
    "balears": [-1, 5, 37.7, 41],
    "med": [-6, 36.33, 30, 46],
}


    # MHWs statistics 

# Basic MHWs statistics
mhws_basic_stats = [
    'count',
    'total_days',
    'duration',
    'total_icum',
    'intensity_max_max',
    'intensity_mean',
]

# Advanced MHWs statistics
mhws_advanced_stats = [
    'intensity_max',
    'intensity_var',
    'rate_onset',
    'rate_decline',
    'intensity_cumulative',
]

# MHWs categories statistics
mhws_categories_stats = [
    'total_days',
    'moderate_days',
    'strong_days',
    'severe_days',
    'extreme_days',
]

# Temperature statistics
mhws_temp_stats = [
    'temp_min',
    'temp_mean',
    'temp_max'
]

# All MHWs statistics
mhws_stats = [
    # Counts of events/days
    'count',
    'total_days',
    'moderate_days',
    'strong_days',
    'severe_days',
    'extreme_days',

    # Cumulative intensity statistics
    'total_icum',
    'intensity_cumulative',

    # Intensity statistics
    'intensity_max_max',
    'intensity_max',
    'intensity_mean',
    'intensity_var',

    # Duration statistic
    'duration',

    # Rate onset/decline statistics
    'rate_onset',
    'rate_decline',

    # Temperature statistics
    'temp_min',
    'temp_mean',
    'temp_max'
]

# Short names for MHWs statistics
mhws_stats_shortname = {
    "count":                "Number of MHW events",
    "total_days":           "Number of MHW days",
    "moderate_days":        "Number of moderate MHW days",
    "strong_days":          "Number of strong MHW days",
    "severe_days":          "Number of severe MHW days",
    "extreme_days":         "Number of extreme MHW days",
    "duration":             "Mean MHW duration",
    "total_icum":           "Total MHW cumulative intensity",
    "intensity_cumulative": "Mean MHW cumulated intensity",
    "intensity_max_max":    "Maximum MHW intensity",
    "intensity_max":        "Mean MHW maximum intensity",
    "intensity_mean":       "Mean MHW intensity",
    "intensity_var":        "Mean MHW intensity variability",
    "rate_onset":           "Mean MHW onset rate",
    "rate_decline":         "Mean MHW decline rate",
    "temp_min":             "Minimum temperature",
    "temp_mean":            "Mean temperature",
    "temp_max":             "Maximal temperature",
}

# Long names for MHWs statistics (from Oliver's code)
mhws_stats_longname = {
    "count":                "Total MHW count per year",
    "total_days":           "Total number of MHW days per year",
    "moderate_days":        "Total number of moderate MHW days per year",
    "strong_days":          "Total number of strong MHW days per year",
    "severe_days":          "Total number of severe MHW days per year",
    "extreme_days":         "Total number of extreme MHW days per year",
    "duration":             "Average MHW duration per year",
    "total_icum":           "Total cumulative intensity over all MHWs per year",
    "intensity_cumulative": "Average MHW \"cumulative intensity\" per year",
    "intensity_max_max":    "Maximum MHW \"maximum (peak) intensity\" per year",
    "intensity_max":        "Average MHW \"maximum (peak) intensity\" per year",
    "intensity_mean":       "Average MHW \"mean intensity\" per year",
    "intensity_var":        "Average MHW \"intensity variability\" per year",
    "rate_onset":           "Average MHW onset rate per year",
    "rate_decline":         "Average MHW decline rate per year",
    "temp_min":             "Minimum temperature per year",
    "temp_mean":            "Mean temperature per year",
    "temp_max":             "Maximum temperature per year",
}

# MHWs statistics units
mhws_stats_units = {
    "count":                "count",
    "total_days":           "days",
    "moderate_days":        "days",
    "strong_days":          "days",
    "severe_days":          "days",
    "extreme_days":         "days",
    "duration":             "days",
    "total_icum":           "°C.days",
    "intensity_cumulative": "°C.days",
    "intensity_max_max":    "°C",
    "intensity_max":        "°C",
    "intensity_mean":       "°C",
    "intensity_var":        "°C",
    "rate_onset":           "°C/days",
    "rate_decline":         "°C/days",
    "temp_min":             "°C",
    "temp_mean":            "°C",
    "temp_max":             "°C",
}

# Colormaps to be used for each MHWs statistics
mhws_stats_cmaps = { # cmo.thermal, Reds, YlOrRd
    "count":                "cmo.ice_r",                 # "Purples",
    "total_days":           "cmo.tempo",                 # "Blues",
    "moderate_days":        "cmo.tempo",                 # "Blues",
    "strong_days":          "cmo.tempo",                 # "Blues",
    "severe_days":          "cmo.tempo",                 # "Blues",
    "extreme_days":         "cmo.tempo",                 # "Blues",
    "duration":             "cmo.matter",                   # "Blues", "cmo.matter"
    "total_icum":           "cmo.solar_r",               # "Oranges",
    "intensity_cumulative": "cmo.solar_r",               # "Oranges",
    "intensity_max_max":    "cmo.amp",                   # "Reds",
    "intensity_max":        "cmo.amp",                   # "Reds",
    "intensity_mean":       "cmo.amp",                   # "Reds",
    "intensity_var":        "cmo.amp",                   # "Reds",
    "rate_onset":           "cmo.speed",                 # "YlOrRd",
    "rate_decline":         "cmo.speed",                 # "YlOrRd",
    "temp_min":             "cmo.thermal",               # "cmo.thermal",
    "temp_mean":            "cmo.thermal",               # "cmo.thermal",
    "temp_max":             "cmo.thermal",               # "cmo.thermal",
}

# Description to add to generated MHWs dataset
mhw_dataset_description = "MHWs yearly statistics computed using the marineHeatWaves " \
        "module for python developped by Eric C. J. Oliver."

# Acknowledgment to add to generated MHWs dataset depending on the original dataset used
cmems_sst_acknowledgment = 'Generated using E.U. Copernicus Marine Service Information, ' \
        'Mediterranean Sea - High Resolution L4 Sea Surface Temperature Reprocessed (DOI: https://doi.org/10.48670/moi-00173)'
cmems_mfc_acknowledgment = 'Generated using E.U. Copernicus Marine Service Information, ' \
        'Mediterranean Sea Physics Reanalysis (DOI: https://doi.org/10.25423/CMCC/MEDSEA_MULTIYEAR_PHY_006_004_E3R1)'

