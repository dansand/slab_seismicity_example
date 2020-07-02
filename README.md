# slab seismicity example


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