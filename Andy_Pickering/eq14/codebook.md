
# Outline of code and processing for EQ14 chi-pod analysis

## Main Scripts

-  _**eq14_patches_paths.m**_ - All paths for data and processed output and figures are set in this script.

-  _**ProcessEq14Cham_AP.m**_ - Process the raw Chameleon profiles and produce mat files with calibrated t,s,p,t' etc. . These files are what the chi-pod method is applied to.

-  _**run_eq14_AP**_ - Run the Chameleon processing to produce 1-m avg profiles of chi and epsilon (using shear probe data). I run this processing again using a smaller fmax because the spectra look like they roll off much lower than the normally assumed 32hz.

- _**Identify_ML_eq14**_ - Identify the mixed layer depth, in order to exclude data where the water column is convectively unstable from further analysis. These data are excluded w/ the functions _**discard_convection_eq14_cham**_ and _**discard_convection_eq14_chi**_.

- _**MakeCasts_eq14.m**_ - Does some pre-processing for chi-pod method, to reduce repeated calculations that take a long time.

-  _**ComputeChi_Chameleon_Eq14**_  - Apply chi-pod method to Chameleon profiles (thermistor data only, not shear probe).

-  _**Make_Overview_Plots_eq14**_ - Make plots for notes (<https://github.com/OceanMixingGroup/Analysis/blob/master/Andy_Pickering/eq14/notes/OverviewNotes.pdf>).




## Params for chi-pod method
When applying the chi-pod method to profiles, a structure **Params** is required, containing the following parameters. The output files are saved in folders named according to these params.
- **fmax** - Max frequency to integrate dT/dt spectrum up to. Determined by where sensor response rolls off, depends on individual thermistor. For eq08, I estimated this to be about 10hz.
- **z_smooth** - The depth interval over which N^2 and dT/dz are smoothed for the chipod calculations.
- **gamma** - Mixing efficiency (Default 0.2).
- **nfft** - # points to use for spectra (default 128)
- **TPthresh** - Threshold for dT/dz spectra power, to avoid doing calculation where signal is very small/noise.
- **resp_corr** - Option to apply response correction to thermistor spectra. 
- **fc** - Cutoff frequency for response correction, if applied.

## Parameters for plots
In Make_Overview_plots_eq14, there are some additional parameters for loading the data to plot.
- **screen_chi** - Discard chi-pod chi and epsilon where log10(epsilon) < -8.5; this is same noise floor used for Chameleon.
- **screen_ml** - Discard data where water column is convectively unstable (**Identify_ML_eq14.m**)
- **Pmin** - Discard all data shallower than Pmin

## Helper Functions

Chameleon processing:
- _**cali_eq14.m**_
- _**tag_file_eq14.m**_

Other 
- _**load_cal_eq14.m**_
- _**gen_mfiles/discard_convection_eq14_cham.m**-
- _**gen_mfiles/discard_convection_eq014_chi.m**-
- _**gen_mfiles/Get_and_bin_profiles.m**_
