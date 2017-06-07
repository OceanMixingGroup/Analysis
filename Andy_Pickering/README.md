# Analyses by Andy Pickering.

This repository contains code and results of analysis done by Andy Pickering as postdoc at OSU. By putting them on github, i'm hoping to:
- Make it easier to share analyses (instead of emailing pdfs etc..)
- Make it easy to update analyses / fix bugs.
- Make it easier for colleagues to see actual code used, giving them a better understanding of the analysis and also spot mistakes/bugs.
- Provide a detailed record and allow others to more easily pick up and continue analyses in future.


## Effect of sampling speed on overturns 
- Motivated by computing overturns from moored profiler (MP) data during IWISE. Wanted to see how the finite sampling speed affected overturn calculations (overturns assume instantaneous profile).
- <https://github.com/andypicke/SimProfiler>

## chi-pod processing routines
- Developed standard Matlab routines for processing CTD-chipod data from cruises. 
- <https://github.com/OceanMixingGroup/mixingsoftware/tree/master/CTD_Chipod>

## EQ14 chameleon/chi-pod analysis 
- Tested chi-pod method by applying it to Chameleon thermistor data, from EQ14 profiles, and comparing to epsilon derived from shear probes.
- Compared CTD-chipod profiles to Chameleon profiles.

## EQ08 chameleon/chi-pod analysis 
- Tested chi-pod method by applying it to Chameleon thermistor data, from EQ08 profiles, and comparing to epsilon derived from shear probes.

## Global patterns from microstructure database
- Looked at relationships between chi and epsilon, estimated mixing efficencies, and epsilon estimated from chi, for several global microstructure datasets.
<https://github.com/OceanMixingGroup/Analysis/blob/master/Andy_Pickering/micro_database/micro_data_notes.pdf>

## Gamma calculations
Calculated mixing efficiency from EQ14 chameleon profiles, using binned data and data computed over patches only. 


