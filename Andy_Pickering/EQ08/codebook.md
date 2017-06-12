
# Outline of code and processing for EQ08 chi-pod analysis

Listed below are the main codes that are run for this analysis, in approximately the order they need to be run in. 

-  _**make_cham_cal_files_eq08_AP.m**_ - Process the raw Chameleon profiles and produce mat files with calibrated t,s,p,t' etc. . These files are what the chi-pod method is applied to.

-  _**run_eq08_avg_AP.m**_ - Run the Chameleon processing to produce 1-m avg profiles of chi and epsilon (using shear probe data). I run this processing again using a smaller fmax because the spectra look like they roll off much lower than the normally assumed 32hz.

- _**Identify_ML_eq08.m**_ - Identify the mixed layer depth, in order to exclude data where the water column is convectively unstable from further analysis. These data are excluded w/ the functions _**discard_convection_eq08_cham**_ and _**discard_convection_eq08_chi**_.

- _**MakeCasts_eq08**_ - Do some pre-processing for chi-pod method to reduce repeated calculations that take a long time.

-  _**ComputeChi_Chameleon_Eq08.m**_  - Apply chi-pod method to Chameleon profiles (thermistor data only, not shear probe).

- _**make_combined_data_files.m**_  - Combine and average chameleon and chi-pod method profiles for different sets of parameters, and save data files that can be loaded when making plots etc.. Loading/combining all the profiles is kind of slow, so good to not have to repeat it when modifying plots or analysis.

-  _**Make_Overview_Plots_eq08.m**_ - Make plots for notes.
