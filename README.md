# ApRES-hydrology

<strong>Matlab code used to produce the results discussed in <i>Surface Meltwater Impounded by Seasonal Englacial Storage in West Greenland</i>. See https://doi.org/10.1029/2018GL079787 for access to manuscript. 
</strong>

<h2>Brief description of filestructure</h2>

<b> Data </b>  
* .mat versions of the data files archived as .csv files in https://doi.org/10.6084/m9.figshare.7104368.v2. 
* Contains indentified internal layer power, bed power, and atmospheric weather station (AWS) temperature and melt data. 

<b> Ice Attenuation Model </b> 
* Matlab scripts used to estimate the expected change in radar attenuation through the top 50m of ice due to seasonal surface temperature variations. 
* <em>tempPowerLoss.m</em> is the wrapper script that uses the iceTempModel and iceAttenuationModel functions to estimate radar attenuation.
  
<b> Water Layer Scattering Model </b>

  Matlab scripts used to both
1. Forward model radar attenuation from observed surface melting for a given water conductivity, englacial layer porosity and pore size.
2. Estimate water layer thickness from observed observed internal layer radar attenuation for a given conductivity, porosity, and pore size

* <em> computeStorageScatt.m </em> is the wrapper script that calculates the estimated water storage from observed internal layer attenuation for a given conductivity, porosity and pore size 
* <em> modelEnglacialAttenuation.m </em> is the wrapper script that forward models radar attenuation from observed AWS melt data for a given conductivity and porosity and a range of pore sizes. This script will loop over the range of pore sizes and estimate the best fit pore size given the radar observed average internal layer attenuation. This script replicates Figure 3 of the manuscript described above
* <em> plotRadarAttenuation.m </em> is the script that plots the observed internal layer power after applying the offset correction between deployment 1 and deployment 2

<h2> Installation Instructions </h2>

1. Fork or download repository
2. Visit https://omlc.org/software/mie/ 
3. Download "Maetzler's MATLAB code for Mie theory" 
4. Place Mie scattering code in current working directory or add to your MATLAB's path


<em> Mie scattering code required to estimate water storage and forward model radar attenuation </em>

<b> Optical Constant Data for changing frequency in volume scattering code <b>
  
  https://mospace.umsystem.edu/xmlui/handle/10355/11599
  https://atmos.washington.edu/ice_optical_constants/
  
 See also: https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2007JD009744

<a href="https://zenodo.org/badge/latestdoi/142612630"><img src="https://zenodo.org/badge/142612630.svg" alt="DOI"></a>
