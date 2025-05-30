########################################################################################################################
########################################################################################################################

###### USER NOTES:
# Tell what the code is about

# This script is runnable

########################################################################################################################
##################################### IMPORTS ##########################################################################
########################################################################################################################

# Basic libraries
import sys
import importlib

########################################################################################################################
##################################### CODES ############################################################################
########################################################################################################################

# All the packages used there, followed by the version used during development
packages = [
    ("matplotlib", "3.10.1"),
    ("numpy", "1.26.4"),
    ("xarray", "2024.2.0"),
    # ("dask", "2024.2.0"),   # Useful for chunkloading with xarray
    ("cartopy", "0.24.1"),
    ("cmocean", "v3.0.3"),
]

# Some printing to decorate the prompt
print()
print(" - Python versions - ")
print("This script is meant to ensure that you have the right versions for Python interpreter and packages.")
print("Versions in parenthesis like (~1.26.4) are the versions used during developement for reference.")
print()
print(f"Python version: {'.'.join(map(str, sys.version_info[:3]))} (~3.10.12)")
print()
print("Packages versions with current environment:")


# Print each package versions
for package_name, dev_version in packages:
    try:
        # If the package is installed, print its version
        module = importlib.import_module(package_name)
        print(f"{package_name}: {module.__version__} (~{dev_version})")
    
    except:
        # If the package is not installed, print it
        print(f"{package_name}: not installed (~{dev_version})")


# Help the user by providing the command to install all the packages
print()
print("Command to install all the packages:")
print(f"pip install {' '.join(package_name for package_name, _ in packages)}")
