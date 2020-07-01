# Readme file for DeltaClassification model

## Requirements

To run these codes, you will need the following software:
* Python 2.7 or earlier (not compatible with Python 3)

The following Python packages are also required:
* matplotlib
* scipy
* numpy
* cPickle
* osgeo
* fiona
* shapely
* utilities
* sklearn
* seaborn
* clusterpy
* itertools
* pandas
* pysal
* collections

## What input is required?
To run this code, the following shape files are required:
* network shapefile, containing the river network extracted from satellite imagery
* island shapefile, containing the land masses or islands of the delta
* patch shapefile, containing the outline of channels

## What does the code do?
The file all.ipynb contains codes run the analysis. From start to finish, the Jupyter Notebook contains code blocks that:
* loads in the shapefiles
* calculate the parameters for the network that both surround and drain the islands
* calculate the base metrics (e.g. perimeter, area, solidity, aspect ratio...)
* calculates maximum distance from the island center to the nearest water body
* estimates minimum, average and maximum widths of all network channels
* evaluates the fractal dimension of each delta island
* creates shapefiles based on the metrics calculated earlier in the code
* saves all metrics to an output file
* generates PCA and GeoSOM results from the island and channel metrics
* plots the U-matrix and dendrogram based on the GeoSOM results
