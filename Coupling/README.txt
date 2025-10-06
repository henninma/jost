*********
Model code, datasets and model output required for coupling (downscaling) of surface mass balance to an evolving ice geometry (dynamics).
The downscaling process is done within the python scripts downscale_smb_spline.py and downscale_smb_spline_future.py, located in the Downscaling directory.
These downscaling scripts are in turn called from issmcoupler.m
issmcoupler.m is in turn called from runme.m
*********

* dem1966_JOB_ref_SMB_100m_UTM32.nc - DTM 1966 at 100 m resolution of ice-surface and land topography, derived from topographic maps and aerial photography.
For details, see the underlying paper Åkesson et al.

* dem2020_JOB_ref_SMB_100m_UTM32.nc - DTM 2020 at 100 m resolution of ice-surface and land topography, derived from https://hoydedata.no/.

* dem_JOB_ref_SMB_1km_UTM32.nc - surface mass balance at 1 km resolution simulated for the present-day (2020) Jostedalsbreen ice cap.

* gl_frac_dummy_JOB_ref_SMB_1km_UTM32.nc - dummy file for surface mass balance at 1 km resolution, used in downscaling process.

* gl_frac_dummy_JOB_ref_SMB_100m_UTM32.nc - dummy file for surface mass balance at 100 m resolution, used in downscaling process.

* monthly_refmb_JOB_1km_1960_2020_UTM32.nc - simulated monthly surface mass balance at 1 km resolution 1960-2020.
Precipitation factor corrected for both ice-covered and ice-free areas. For details, see Sjursen et al. 2025 (doi: https://doi.org/10.1017/aog.2024.41) 

* monthly_refmb_JOB_1km_1960_2020_UTM32_Jan24_Pcorr_glacieronly -  simulated monthly surface mass balance at 1 km resolution 1960-2020.
Precipitation factor corrected for ice-covered areas. For details, see Sjursen et al. 2025 (doi: https://doi.org/10.1017/aog.2024.41) 

* monthly_refmb_JOB_1km_2021_2100_UTM32.nc - simulated monthly surface mass balance at 1 km resolution 2021-2100,
using temperature and precipitation fields from the climate model combination ECEARTH/CCLM.
For more details, see the underlying paper Åkesson et al.
