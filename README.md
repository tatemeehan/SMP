# SMP

_processSMP.m_ is the main script. 

The process is parallel enabled, which will speed up the inversion loop when the control flag isParallel = 1. The default number of threads is 4, as this is fairly common on a laptop workstation.

After running this script, the user will be prompted to select the directory containing the many SMP profiles desired to be processed. Alternatively, the user may skip this step be seting the control flag isUI = 0 and inputing the data directory where instructed.

## Processing Workflow
Processed SMP data is stored in a cell array containing the SMP data structure of each file called _SMP_. Inversion results are stored in the cell structure _invSMP_.
### Quality Control
The data directory will be searched for the .PNT or .pnt binary SMP files & a .csv file containing the Quality Control metadata. If the QC.csv file is not found, the user will be prompted to QC the data. The SMP profiles will be loaded and plotted in a docked matlab window. The user will then be asked to identify the SMP profile as "Good" - meaning there is no error in the data - "Drift" - meaning an obvious linear drift exists in the profile - "Dry Run" - meaning this file was used as air calibration and will be discarded - "Bad" - an obvious, fatal flaw has occured with the profile and it will be discarded from the process. After completing the QC process, the .csv file will be created and this step will not need to occur again for the given data.
### Automatic Snow Surface Detection
The snow surface is a spike (often subtle) in the force profile following the initial 10 cm or so where the penetrometer is driving through the air. The detection algorithm was designed as such. The variance of the SMP force profile is calculated in a moving window 1 mm in length, the variance is then smoothed by moving average window of 0.25 mm. The initial 5 mm of the force profile are then averaged to create a variance threshold for the air signal (used a proxy for the instrument noise level). The snow surface is reliably identified as the first index which exceeds the threshold by 7 standard deviations.
### Drift Corection
If the signal has an appreciable drift it will be corrected by a linear approximation. The SMP signal above the snow surface is fit with a linear function by least-squares regression. This linear segment is then extrapolated to the entire length of the SMP profile. The drift correction is applied by subtracting this linear drift function from the SMP penetration force profile.
### Depth Correction
The force data above the snow surface is removed and the depth axis is reconfigured to represent the shortened force profile - where zero is the snow surface.
### Inversion for Micromechanical Snow Properties
The inversion strategy of _Marshall & Johnson (2009) doi:10.1029/2009JF001269_ is implemented in this workflow. The microstructural parameters, the rupture force $f$, the deflection at rupture $\delta$, and the characteristic length $l$ are first estimated. Micromechanical properties of the snow are then derived from these microsctuctural building blocks.
