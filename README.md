# PRJ_Error_eeg
Preprocessing, analysis, and modeling of behavior and EEG data for sequential prediciton error EEG Target Time study.
The goal of the study is to understand the relationship between different learning-related computations
are represented in the evoked potentials recorded in scalp electroencephalography (EEG).
Dataset includes 32 good EEG datasts (41 total collected).

Manuscript is currently under review (9/3/20); written by Colin W. Hoy.

## Dependencies
OS: MacBook Pro running OS 10.13.6; MATLAB version R2017b; Python 2.7 (not tested on any other platforms)
  - External toolboxes:
    - Fieldtrip: <http://www.fieldtriptoolbox.org/>
    - CircStat Toolbox: <https://github.com/circstat/circstat-matlab>
    - FMA Toolbox: <https://github.com/michael-zugaro/FMAToolbox>
  - see colin_PRJ_Error_py2.7.yml for python environment dependencies, which can be replicated in a Conda (or similar) environment (see <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file>)

## Overview of Code
Code execution to reproduce all analyses in the paper can be found in run_lists.
These run scripts contain parameters used to call each function. All analyses should run in under an hour, depending on computing capacity.
1. run00_behavior_preprocessing.m
    - BHV scripts preprocess and analyze behavioral log files.
    - SBJ00, SBJ01, and SBJ02 scripts preprocess EEG data.
2. run01_ERP_scripts.m
    - SBJ03 scripts preprocess and plot event-related potentials (ERPs).
3. run02_TFR_scripts.m
    - SBJ05 scripts preprocess, plot, and model time-frequency representations (TFRs) of EEG data.
4. run03_RL_model_results.m
    - SBJ04 scripts model behavioral and ERP data and plot results.
    - SBJ05 scripts model TFR (power and phase) data and plot results.
5. run04_RL_mdoel_point_estimates.m
    - SBJ04 scripts to compute, model, and plot point estimates of FRN (mean window and peak-to-peak analyses)

## Execution Parameters
Specific parameters for different preprocessing and analysis scripts are loaded as options.
Generally, these options are written as executable MATLAB code that is specified when calling
and analysis script, which then run the relevant code to load the parameters inside that script.

- SBJ_lists: text lists of groups of subjects
  - 'good1' was the initial cohort used to develop, test, and finalize all analyses and parameters
  - 'good2' was the replication cohort
  - **'goodall'** includes both cohorts
    - This was used in the paper to report the most representative results.
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

## Example Data
All raw datasets will be provided upon publication. I originally attempted to provide two example datasets in this repository, but the file sizes were too large and could not be merged. Anyone looking to use the data before publication is welcome to request it from hoycw (at) berkeley (dot) edu.
