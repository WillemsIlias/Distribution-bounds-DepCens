# Distribution-bounds-DepCens
This directory contains all the computer code used in "Bounds for the regression parameters in dependently censored survival models". The directory is organized as follows:

* The implementation of the methodology in R can be found in the following files:
	* chronometer.R:			Used to obtain run times of simulations.
	* EstimationAlgorithmBei.R:		Implementation of the test by Bei (see manuscript for reference).
	* lowLevelFunctions.R:			Main estimation function and its various components.
	* searchStrategies.R:			Functions related to root finding algorithms.
	* simulationAnalysisFunctions.R:	Functions used in analyzing the simulation results.
	* simulationFunctions.R:		Functions used to run all simulations.
	* Data application functions.R:		Functions used to run the data applications.

* All simulations have been carried out making use of the computing clusters in the Flemish Supercomputer Center (VSC), and hence the code to carry out the simulations is written accordingly. All files containing the computer code for executing simulations are named with the prefix "simMain" or "simMain_add", possibly followed by a specification of which simulation it pertains to (e.g. "MoreIF", alluding to the simulation using more instrumental functions). Lastly, a suffix "_wrapper" is appended if the script is used in executing the simulation, while a suffix "_analysis" indicated that the script is used in analyzing the results obtained from the corresponding "_wrapper" script.

* Each "_wrapper" script executes a simulation according to a specified simulation setting. The settings with which to run each of the "_wrapper" files are stored in .csv-files and are organized in folders with names starting by "simMainDesigns".

* The results of the simulations are stored in the folders prefixed with "Simulations_". These folders are accessed by the scripts suffixed with "_analysis" (see also the second point). We provide them here in .zip format, and hence they should be unzipped before they can be used. [Note: Due to exceedance of the maximum file size, these files could not be uploaded.]

* The script executing the data analysis on pancreatic cancer is called "DataApplicationPancreasWrapper_new.R". The settings with which it was run ("DataApplicationPANCREAS_settings.csv"), the preprocessed data ("Pancreas_preprocessed.csv") and the results of the analyses ("PANCREAS results new") are collected in the folder "Data application > PANCREAS".

* The script executing the NLSY data analysis is called "DataApplicationNLSYWrapper.R". The settings with which it was run ("criteria2check_NLSY.csv"), the preprocessed data ("NLSY_unemp_spells.csv") and the results of the analyses ("Results") are collected in the folder "Data application > NLSY79".

Note that running the simulations or the data applications locally might take a lot of time. The file "Fail codes.txt" documents the meaning of possible failure codes returned by the main estimation algorithm ("pi.surv").


