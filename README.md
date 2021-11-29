# slab seismicity example

## Notes

As of Sep 2020, this example is working with py3, using a conda environment (py3_mapping). 

There is now a master branch (python2) and a branch called py3

the python scripts were moved out of this directory to `~/projects/shared_tools`

the py3 branch will look for them at that location


## To do

Replace gdal raster analysis with rasterio (i.e. intepolation). 

Simplify/improve the `create_profile_lines` functionality. Try to replace pygplates with pyproj.Geod, which has great circle tools.


## overview.ipynb

* Overview of the analysis
* python requirements
* notes

## create_profile_lines.ipynb

* define the domain
* create profile lines
* interpolate slab data along the profile lines
* save required data to `./output_data`

## slab_seismicity.ipynb

* read in earthquake catalogs using obspy
* find common events in gCMT amd ISC-EHB catalogs
* calculate earthquake location in the 2d planes defined by profile lines
* calculate earthquake locations in slab reference sytem
* demonstrates plotting of rotated principal axes