![version](https://img.shields.io/badge/version-0.1.0-blue)
[![HitCount](http://hits.dwyl.com/{AlexandreWadoux}/{soilspec}.svg)](http://hits.dwyl.com/{AlexandreWadoux}/{soilspec})
[![Project Status: Active The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![minimal R version](https://img.shields.io/badge/R%3E%3D-3.2.4-6666ff.svg)
[![CRANstatus](https://www.r-pkg.org/badges/version/soilspec)](https://cran.r-project.org/package=soilspec)

<a href="url"><img src="https://github.com/AlexandreWadoux/soilspec/blob/master/drawing.png" align="centre" height="450" width="900" ></a>

# Soil Spectral Inference with R

Data and functions for the book *Soil Spectral Inference with R* - *Analysing Digital Soil Spectra Using the R Programming Environment*. 

## Short description
### Functions

* `chBLCext`: function to fit a convex hull to the region of interest. The regions of interest are absorbance or reflectance of some secondary clay minerals or iron oxides from vis-NIR spectra. The function contains three arguments and returns five values. The function argument description is obtained by `?soilspec::chBLCext`.

* `eval`: function to compute of a set of accuracy measures between observed and predicted continuous values or classes. The user must specify the type of variable, either quantitative or categorical. For both types of variable, a set of accuracy measures is reported. The list of accuracy measures is provided by writing `?soilspec::eval`. The function contains three arguments and returns a number of values equal to the number of accuracy measures computed. 

* `myImagePlot`: function to plot an image from a matrix. This simple function takes as single argument a matrix and returns a plot of this matrix. The functions details are accessed by writing `?soilspec::myImagePlot`.

* `spectra2colour`: this function converts spectra reflectance into RGB and Munsell colours. The functions take a single argument as input, a matrix or data frame of the spectra, and returns the colour for both RBG and Munsell charts. The function returns a data frame where each row is the spectrum followed by the colour. The function description is obtained by `?soilspec::spectra2colour`.

* `css`: function to determine the optimal number of spectra to be sent to the laboratory for soil analysis. This function works by comparing the probability density function (pdf) of the population to that of the sample set to assess the representativeness of the sample. The two pdfs are compared based on the mean Euclidean distance (msd). The functions contains eight arguments and returns three values. The function argument description is obtained by `?soilspec::css`. 

### Datasets

* `datEPO`: subset of the Geeves soil survey data in which the soil samples were scanned with different level of water content. The dataset is a list with four data frames. Description of each data frame is provided in the R documentation by writting `?soilspec::datEPO`.  

* `datsoilspc`: spectra and associated values of laboratory analysed soil properties from the soil samples described in Geeves et al., (1994). The data provided is a data frame with five columns and 391 rows. The first four columns contain values of the clay, silt, sand and total carbon. The fifth columns is a matrix with the infrared spectra. The documentation can be accessed by writting `?soilspec::datsoilspc`. 

* `datStand`: vis-NIR spectra of a standard material as described in Ben-Dor et al., (2015). The data are provided as a list with four data frames. The documentation and description of each of them is provided by writting `?soilspec::datStand`. 

* `mineralRef`: a data frame containing specific mineral compounds measured by laboratory spectrometers. The data frame contained 13 columns. The first is the wavelength and the others are 12 spectra from mineral compounds. The full description is provided by writting `?soilspec::mineralRef`. 

* `rutherglenNIR`: Small set of vis-NIR spectra collected from a soil profile. The documentation is provided by writting `?soilspec::rutherglenNIR`. 

* `specSoilCol`: Small set of vis-NIR spectra with value of the soil horizon and soil type. In addition to the spectra, pictures of the soil samples were made and provided in the book. This enables the comparison of the colour retrieved from the spectra and that from the picture of the soil samples. The documentation is provided in `?soilspec::specSoilCol`. 

## Installation
The `soilspec` package is available in this GitHub repository only. The package is likely to be updated in the future. It is recommended to re-install the package is you didn't use it for a long time. If any issue, send an email to [Alexandre Wadoux](mailto:alexandre.wadoux@sydney.edu.au) or open an [issue](https://github.com/AlexandreWadoux/soilspec/issues). 

The package is installed using the following command. 
```R
# check if the devtools package is already installed
if (!require("devtools")) install.packages("devtools")

# install the soilspec package from GitHub
devtools::install_github("AlexandreWadoux/soilspec")
```

## Help
Help are obtained using the R documentation provided in the package `?soilspec::<function>`. 
The author can be contacted by email: alexandre.wadoux@sydney.edu.au. 

## Reference
Wadoux et al., ---
