
# Outline of code and processing for EQ08 chi-pod analysis

-  _**make_cham_cal_files_eq08_AP**_ - Process the raw Chameleon profiles and produce mat files with calibrated t,s,p,t' etc. . These files are what the chi-pod method is applied to.

-  _**run_eq08_avg_AP**_ - Run the Chameleon processing to produce 1-m avg profiles of chi and epsilon (using shear probe data). I run this processing again using a smaller fmax because the spectra look like they roll off much lower than the normally assumed 32hz.

- _**Identify_ML_eq08**_ - Identify the mixed layer depth, in order to exclude data where the water column is convectively unstable from further analysis.

-  _**ComputeChi_Chameleon_Eq08**_  - Apply chi-pod method to Chameleon profiles (thermistor data only, not shear probe).

-  _**Make_Overview_Plots_eq08**_ - Make plots for notes.
