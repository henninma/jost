# Åkesson et al. 2025: supporting code
Model code supporting the paper about Jostedalsbreen ice cap by Åkesson et al. 2025

Metadata and software details for model code accompanying the publication
Åkesson, H., Sjursen, K. H., Vikhamar Schuler, T., Dunse, T., Andreassen, L. M., Kusk Gillespie, M., Robson, B. A., Schellenberger, T., and Clement Yde, J.: Recent history and future demise of Jostedalsbreen, the largest ice cap in mainland Europe, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2025-467, 2025. REFERENCE TO BE UPDATED

Please contact the corresponding author henning.akesson [at] geo.su.se for any scientific or software-related questions.

1. Overview
This repository contains code used to prepare and launch simulations of Jostedalsbreen ice cap in Norway, using the 'Ice-sheet and Sea-level System Model' (ISSM). The model code needed to run such simulations is freely available at https://issm.jpl.nasa.gov/download/.

3. Matlab requirements
The scripts have been run on Matlab R2024a and R2024b, but should work on other versions as well.

5. Other requirements
To be able to launch the scripts included here, the user must first download and compile the model code of the 'Ice-sheet and Sea-level System Model' (ISSM) at https://issm.jpl.nasa.gov/download/.
The necessary underlying datasets for ice geometry, ice velocities and surface mass balance are listed in the 'Code and data availability' section in the published paper [DOI TO BE ADDED]. 


7. Description of files included

- runme.m - 
Main script to prepare and launch simulations. 
Sets up the model mesh, necessary parameterizations of model processes
and climate forcing, such as basal motion, surface mass balance, and launches and saves 
the model simulations.

- Scripts/Jost.par - 
Underlying script, used by runme.m, containing parameterizations of some 
key processes as well as some model settings required.

- .exp files -
Contain coordinates delineating the model domain, as well as some areas
where specific constraints or boundary conditions are set within runme.m.
Coordinates are given in UTM format. The projection is UTM32N.
The .exp files have the format

## Name:Jost
## Icon:0
# Points Count Value
84 1.000000
# X pos Y pos
383055.7773288192 6843732.3558886303
386674.9498249681 6847544.2450076826
388192.2207827155 6848420.7355066249
