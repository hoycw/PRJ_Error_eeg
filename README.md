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
  - see colin_PRJ_Error_py2.7.yml for python einvironment dependencies

## Overview of Code
Code execution to reproduce all analyses in the paper can be found in run_lists.
These run scripts contain parameters used to call each function.
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
Two sample datasets are provided in the folder demo_data. These two subjects were chosen as representative of good (EEG01) and medium (EEG06) quality data. Upon publication, all datasets will be made available in a public repository.

### Preprocessing Demo Data
Raw data files are provided in demo_data/EEG**/00_raw/.
To preprocess these datasets, create new SBJ_vars for them by copying the originals and modifying the file paths. Please note that re-running preprocessing scripts will require changing specific ICA component and trial rejection indices.

### Analyzing Demo Data
Preprocessed EEG data files (using proc_id = 'eeg_full_ft') are provided in demo_data/EEG**/02_preproc/.
Corresponding behavioral data are provided in demo_data/EEG**/03_events/.
These datasets have been cleaned and can be used to run ERP, TFR, and modeling scripts to reproduce analyses. By using both subjects, group scripts can be tested, in which case the SBJ_lists/demo.sbj should be used as SBJ_id when calling analysis scripts.
Note that path names in analysis and plotting scripts will need to be adjusted depending on directory structure.
