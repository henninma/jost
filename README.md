# Åkesson et al. 2025, supporting code and assets
Model code and metadata supporting the paper about Jostedalsbreen ice cap by Åkesson et al. 2025 in The Cryosphere.

Åkesson, H., Sjursen, K. H., Vikhamar Schuler, T., Dunse, T., Andreassen, L. M., Kusk Gillespie, M., Robson, B. A., Schellenberger, T., and Clement Yde, J.: Recent history and future demise of Jostedalsbreen, the largest ice cap in mainland Europe, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2025-467, 2025. REFERENCE TO BE UPDATED

Please contact the corresponding author Henning Åkesson (henning.akesson [at] geo.uio.no) for any scientific or software-related questions.

## Overview
This repository contains code used to prepare and launch simulations of Jostedalsbreen ice cap in Norway, using the 'Ice-sheet and Sea-level System Model' (ISSM). The model code needed to run such simulations is freely available at https://issm.jpl.nasa.gov/download/.

## Videos of future evolution
The following animations are supplied along with the article, and can be downloaded from the videos folder (permalink Zenodo TO BE ADDED). These are available in .avi and .mp4 formats. The latter works best with PowerPoint.

– Supplementary video 1: Modelled ice-thickness evolution 2021–2100 (ECEARTH/CCLM, RCP4.5).

– Supplementary video 2: Modelled ice-thickness evolution 2021–2100 (ECEARTH/CCLM, RCP8.5).

– Supplementary video 3: Modelled ice-thickness evolution 2101–2300 (Commit4.5).

– Supplementary video 4: Modelled ice-thickness evolution 2101–2300 (Commit8.5). 

## Software requirements
The scripts have been run on Matlab R2024a and R2024b, but should work on other versions as well. If you'd like to run coupled simulations of surface mass balance and ice dynamics, you will also need a Python virtual environment, see instructions below.

## Other requirements
To be able to launch the scripts included here, the user must first clone and compile the model code of the 'Ice-sheet and Sea-level System Model' (ISSM) at https://issm.jpl.nasa.gov/download/.
The necessary underlying datasets for ice geometry, ice velocities and surface mass balance are listed in the 'Code and data availability' section in the published paper [DOI TO BE ADDED]. If additional data is needed to support your simulations, please contact the corrsonding author Henning Åkesson (henning.akesson [at] geo.uio.no).

## Description of files included
- runme.m - 
Main script to prepare and launch simulations. 
Sets up the model mesh, necessary parameterizations of model processes
and climate forcing, such as basal motion, surface mass balance, and launches and saves 
the model simulations.

- Scripts/Jost.par - 
Underlying script, used by runme.m, containing parameterizations of some 
key processes as well as some model settings required.

- Coupling/createNetcdfSMB.m -
Handles writing of surface mass balance to netcdf, needed during the coupling process in simulations.

- Coupling/issmcoupler.m -
Handles the coupling of surface mass balance to ice dynamics (evolving geometry).

- Coupling/Downscaling/downscale_smb_spline_future.py -
Handles downscaling of reference surface mass balance to the current ice geometry during the coupling in future simulations.

- Coupling/Downscaling/downscale_smb_spline.py -
Handles downscaling of reference surface mass balance to the current ice geometry during the coupling in historical simulations.

- Coupling/Downscaling/base_files/ -
Datasets needed during the downscaling of surface mass balance. For details, see readme file within this directory. 

- .exp files -
Contain coordinates delineating the model domain, as well as some areas
where specific constraints or boundary conditions are set within runme.m.
Coordinates are given in UTM format. The projection is UTM32N.
The .exp files have the format

```## Name:Jost
## Name:Domain
## Icon:0
# Points Count Value
84 1.000000
# X pos Y pos
383055.7773288192 6843732.3558886303
386674.9498249681 6847544.2450076826
388192.2207827155 6848420.7355066249
```


## Coupled simulations of surface mass balance and ice dynamics (Python + Matlab)
To run coupled simulations of surface mass balance and ice dynamics (evolving geometry) you will, in addition to Matlab, need a Python virtual environment, which can be used with Matlab's Python Interface. 

The user can choose how often surface mass balance and ice dynamics should be coupled. This is set using the variable ```coupling_interval``` in the beginning of the script runme.m.

The downscaling process of surface mass balance to the evolving ice geometry is done within the python scripts downscale_smb_spline.py and downscale_smb_spline_future.py. These downscaling scripts are in turn called from issmcoupler.m, which in turn is called from runme.m.

To enable Matlab to call python scripts used to downscale surface mass balance to an evolving ice geometry (dynamics), we'll need to create a Python virtual environment, see [this guide](https://se.mathworks.com/support/search.html/answers/1750425-python-virtual-environments-with-python-interface.html?fq%5B%5D=asset_type_name%3Aanswer&page=1). For example:

Assume your python path is:
```/opt/software/easybuild/software/Python/3.9.6-GCCcore-11.2.0/bin/python```

Make a Python virtual environment in your project directory:

```mimi:jost$ /opt/software/easybuild/software/Python/3.9.6-GCCcore-11.2.0/bin/python -m venv /uio/hypatia/geofag-personlig/geohyd-staff/henninma/issm/trunk-jpl/projects/jost/py39```

Make sure pip is upgraded:

```/uio/kant/geo-geofag-u1/henninma/issm/trunk-jpl/projects/jost/py39/bin/python -m pip install --upgrade pip```

Source environment before installation of python packages:

```source /uio/hypatia/geofag-personlig/geohyd-staff/henninma/issm/trunk-jpl/projects/jost/py39/bin/activate```

Now your virtual environment called ```py39``` should be activated.

Install python packages needed:

```(py39) /home/username$ python -m pip install numpy xarray scipy datetime```

You will also need to install input/output tools for xarray, [see here](https://docs.xarray.dev/en/stable/getting-started-guide/installing.html)

```python -m pip install "xarray[io]"```

For netcdf into using xarray, see also [here](https://docs.xarray.dev/en/stable/generated/xarray.Dataset.to_netcdf.html)

If you run into trouble during the creating of the Python virtual environment for Matlab, please refer to [this guide](https://se.mathworks.com/support/search.html/answers/1750425-python-virtual-environments-with-python-interface.html?fq%5B%5D=asset_type_name%3Aanswer&page=1).

Once you have a working Python virtual environment, you'll need to activate it everytime you would like to run coupled simulations, and **before** you start Matlab. You can add an alias to your .bashrc like this:

```alias py39='source /uio/hypatia/geofag-personlig/geohyd-staff/henninma/issm/trunk-jpl/projects/jost/py39/bin/activate'```

Modify the path to your python virtual environment accordingly. This means you can then just do ```py39``` in your terminal to activate the virtual environment.


