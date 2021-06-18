# PRJ_Error_eeg: Sequential Prediction Errors in Reward Processing EEG
Preprocessing, analysis, and modeling of behavior and EEG data for sequential prediciton error EEG Target Time study.
The goal of the study is to understand the relationship between different learning-related computations
are represented in the evoked potentials recorded in scalp electroencephalography (EEG).
Dataset includes 32 good EEG datasets (41 total collected) and 22 good behavioral ratings datasets (24 total collected).
Note that "TT" refers to the Target Time task, and "OB" refers to the three tone Oddball task.
Key predictors: Expected value (EV), RPE value (sRPE), RPE magnitude (uRPE), probability (Lik)

## Manuscript: 
Manuscript has been accepted for publication at Communications Biology and is in press.
First round of revisions submitted April 16, 2021, and second round of revisions submitted June 5, 2021.
Initial biorxiv preprint: https://www.biorxiv.org/content/10.1101/2020.10.20.347740v1
Code written by Colin W. Hoy.

## Dependencies
OS: MacBook Pro running OS 10.13.6; MATLAB version R2017b; Python 2.7 (not tested on any other platforms)
  - External toolboxes:
    - Fieldtrip: <http://www.fieldtriptoolbox.org/>
    - CircStat Toolbox: <https://github.com/circstat/circstat-matlab>
    - FMA Toolbox: <https://github.com/michael-zugaro/FMAToolbox>
  - see colin_PRJ_Error_py2.7.yml for python environment dependencies, which can be replicated in a Conda (or similar) environment (see <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file>)

## Overview of Code
Code execution to reproduce all analyses in the paper can be found in run_lists.
These run scripts contain parameters used to call each function. All individual analyses should run in under an hour, depending on computing capacity.
1. run00_behavior_preprocessing.m
    - BHV scripts preprocess and analyze behavioral log files.
    - SBJ00, SBJ01, and SBJ02 scripts preprocess EEG data.
2a. run01a_ERP_preprocessing_TT.m
    - SBJ03 scripts preprocess and plot event-related potentials (ERPs).
2b. run01b_ERP_preprocessing_odd.m
    - ODD01, ODD02a, ODD02c versions of preprocessing
    - SBJ03 scripts preprocess and plot event-related potentials (ERPs).
2c. run01c_difference_waves.m
    - SBJ03 scripts compute and plot difference waves.
3a. run02a_TFR_preprocessing_TT.m
    - SBJ05 scripts preprocess, plot, and model time-frequency representations (TFRs) of TT EEG data.
3b. run02b_TFR_preprocessing_OB.m
    - SBJ07 scripts preprocess, plot, and model time-frequency representations (TFRs) of OB EEG data.
4. run03_RL_model_results.m
    - SBJ04 scripts model behavioral and ERP data and plot results.
    - SBJ05 scripts model TFR power data and plot results.
5. run04_RL_model_point_estimates.m
    - SBJ04 scripts to compute, model, and plot point estimates of FRN (mean window and peak-to-peak analyses)
6. run05_OB_TT_feature_comparison.m
    - SBJ06 scripts to extract ERP features from TT and OB, run feature correlations
7. run06_rating_comparisons.m
    - preprocessing and analysis of behavioral rating data

## Execution Parameters
Specific parameters for different preprocessing and analysis scripts are loaded as options.
Generally, these options are written as executable MATLAB code that is specified when calling
and analysis script, which then run the relevant code to load the parameters inside that script.

- SBJ_lists: text lists of groups of subjects
  - 'good1' was the initial cohort used to develop, test, and finalize all analyses and parameters
  - 'good2' was the replication cohort
  - **'goodall'** includes both cohorts
    - This was used in the paper to report the most representative results.
  - 'goodOB' was the subjects with Oddball data passing QA procedures
  - 'goodEEG*' are subject lists with both TT and OB data but without OB QA
  - 'ratings_good' was the subjects with behavioral rating data passing QA procedures
- SBJ_vars: contains subject specific information
  - Raw data file names
  - EEG information (e.g., channel labels)
  - Preprocessing information (e.g., bad channels, trials, independent components, etc.)
- questionnaire_vars: manual coding of survey data for each participant
- proc_vars: preprocessing parameters for artifact rejection
- an_vars: analysis parameters for ERP and TFR computations
  - Channel selection, epoching, filtering, baseline correction, etc.
- stat_vars: statistics parameters for modeling behavioral, ERP, and TFR data
  - Model regressors, trials/conditions, regression style, epoching and averaging metrics, etc.
- plt_vars: plotting parameters
  - Epochs, styles, markers, significance, legends, etc.
- feat_vars: analysis parameters to extract ERP features from Oddball and Target Time data for correlations
  - measure of ERP activity, conditions, channels, epochs, peak selection, etc.

## Example Data
All raw datasets will be provided upon publication. I originally attempted to provide two example datasets in this repository, but the file sizes were too large and could not be merged. Anyone looking to use the data before publication is welcome to request it from hoycw (at) berkeley (dot) edu.
