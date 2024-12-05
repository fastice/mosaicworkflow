# mosaicWorkflow

This repository contains the Python code used to produce the velocity and image mosaics for the Greenland Ice Mapping Project (GrIMP). Mosaics are produced by call a set of C-languge programs, which can be accessed from [here](https://github.com/fastice/mosaicSource).

This is production code. Most users will be better served by an extensive set of [Jupyter Notebooks and Python packages](https://github.com/fastice/GrIMPTools) to work with the GrIMP [products](https://nsidc.org/data/measures/grimp).

## Mosaicking Program Descriptions

The following command line applications generate velocity and image mosaics by calling the C-language command-line [modules](https://github.com/fastice/mosaicSource). Each of the C-routines is single-threaded. To apply parallism, the Python modules chunk up so that it can be processed with multiple threads, which are then reassembled to produce the final product. 

Outputs are all the final archive-ready products that get delivered to NSIDC.

**setupquarters**: This program produces single cycle, monthly, quarterly, or annual mosaics from a combination of SAR and Optical data. It is often not called directly but is instead called from **makemosaic**

**makemosaic**: Performs preprocessing and post processing task to create a time series of mosaics using **setupquarters**. It performs functions like selecting the appropriate mask for the time period, creating an input file that contains only data relevant to the specified time range, and ensuring the date ranges conform to the specified product ranges.

**setupimagemosaic**: Produces a single image mosaic or a mosaic with sigma0 and gamma0 bands. As with the velocity mosaics, the product is chunked up, farmed to several threads, and then reassembled.

**makeimagemosics** Peforms the necessary preprocessing steps to produce a series of image mosaics using **setupimagemosaic**.

## Other Utilities

**simoffsets** Simulate range and azimuth offsets using a velocity map and DEM.

**ROFFtoGrimp** Pre-process a NISAR *ROFF* product to remove outliers and then merge the three offset bands into a single offset product as a set of binary files.

**RUNWtoGrimp** Pre-process a NISAR *RUNW* product to apply masks (e.g., connected components) and corrections (e.g., ionosphere) and output the result as binary file along with a geojson file containing relevant metadata. 

## For Further Information

Please address questions to ![](https://github.com/fastice/GrIMPTools/blob/main/Email.png).
