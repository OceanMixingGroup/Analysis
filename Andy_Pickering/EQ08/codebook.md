
# Outline of code and processing for EQ08 chi-pod analysis

## Main Scripts
Listed below are the main codes that are run for this analysis, in approximately the order they need to be run in. 

-  _**eq08_patches_paths.m**_ - (in gen_mfiles/) Sets all the paths for data and output.

-  _**make_cham_cal_files_eq08_AP.m**_ - Process the raw Chameleon profiles and produce mat files with calibrated t,s,p,t' etc. . These files are what the chi-pod method is applied to.

-  _**run_eq08_avg_AP.m**_ - Run the Chameleon processing to produce 1-m avg profiles of chi and epsilon (using shear probe data). I run this processing again using a smaller fmax because the spectra look like they roll off much lower than the normally assumed 32hz.

- _**Identify_ML_eq08.m**_ - Identify the mixed layer depth, in order to exclude data where the water column is convectively unstable from further analysis. These data are excluded w/ the functions _**discard_convection_eq08_cham**_ and _**discard_convection_eq08_chi**_.

- _**MakeCasts_eq08**_ - Do some pre-processing for chi-pod method to reduce repeated calculations that take a long time.

-  _**ComputeChi_Chameleon_Eq08.m**_  - Apply chi-pod method to Chameleon profiles (thermistor data only, not shear probe).

- _**make_combined_data_files.m**_  - Combine and average chameleon and chi-pod method profiles for different sets of parameters, and save data files that can be loaded when making plots etc.. Loading/combining all the profiles is kind of slow, so good to not have to repeat it when modifying plots or analysis.

-  _**Make_Overview_Plots_eq08.m**_ - Make plots for notes.


## Params for chi-pod method
When applying the chi-pod method to profiles, a structure **Params** is required, containing the following parameters. The output files are saved in folders named according to these params.
- **fmax** - Max frequency to integrate dT/dt spectrum up to. Determined by where sensor response rolls off, depends on individual thermistor. For eq08, I estimated this to be about 10hz.
- **z_smooth** - The depth interval over which N^2 and dT/dz are smoothed for the chipod calculations.
- **gamma** - Mixing efficiency (Default 0.2).
- **nfft**
- **TPthresh**
- **resp_corr** - Option to apply response correction to thermistor spectra. 
- **fc** - Cutoff frequency for response correction, if applied.

## Parameters for plots
In Make_Overview_plots_eq08, there are some additional parameters for loading the data to plot.
- ** screen_chi** - Discard chi-pod chi and epsilon where log10(epsilon) < -8.5; this is same noise floor used for Chameleon.
- **screen_ml** - Discard data where water column is convectively unstable (**Identify_ML_eq08.m**)
- **Pmin** - Discard all data shallower than Pmin

## Helper Functions

- _**cali_eq08.m**_
- _**tag_file_eq08.m**_
- _**modify_header_eq08.m**_
- _**load_cal_eq08.m**_
- _**discard_convection_eq08_cham.m**-
- _**discard_convection_eq08_chi.m**-
