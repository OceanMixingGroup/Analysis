
# Outline of code and processing for EQ14 chi-pod analysis

-  _**eq14_patches_paths.m**_ - All paths for data and processed output and figures are set in this script.

-  _**ProcessEq14Cham_AP.m**_ - Process the raw Chameleon profiles and produce mat files with calibrated t,s,p,t' etc. . These files are what the chi-pod method is applied to.

-  _**?**_ - Run the Chameleon processing to produce 1-m avg profiles of chi and epsilon (using shear probe data). I run this processing again using a smaller fmax because the spectra look like they roll off much lower than the normally assumed 32hz.

- _**Identify_ML_eq14**_ - Identify the mixed layer depth, in order to exclude data where the water column is convectively unstable from further analysis. These data are excluded w/ the functions _**discard_convection_eq14_cham**_ and _**discard_convection_eq14_chi**_.

- _**MakeCasts_eq14.m**_ - Does some pre-processing for chi-pod method, to reduce repeated calculations that take a long time.

-  _**ComputeChi_Chameleon_Eq14**_  - Apply chi-pod method to Chameleon profiles (thermistor data only, not shear probe).

-  _**Make_Overview_Plots_eq14**_ - Make plots for notes.
