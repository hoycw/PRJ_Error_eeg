#!/bin/bash
# Colin ran this script in his local directories on Jan 16 2020

# Sheila comments when writing this script (Jan 14 2020):
# List of Commands to Run:
#   **I compiled a list of things I deleted from the cluster but 
#   unfortunately did not have the foresight to include the specific sub-directories 
#   they used to reside in, so my apologies if the cd commands error at any point!

cd /data/
# Old Files - all SBJ have
rm */03_events/*_behav_eeg_full_ft_clean.mat
rm */03_events/*_behav02a_eeg_full_ft_clean.mat
rm */03_events/*_behav02b_eeg_full_ft_clean.mat

# Old Files - all EEG SBJ have
rm */03_events/*_behav_odd_full_ft_clean.mat
rm */03_events/*_behav02a_odd_full_ft_clean.mat
rm */03_events/*_behav02b_odd_full_ft_clean.mat

# Old Files - only some SBJs have
rm */03_events/*_behav_oddball_eeg_full_ft_clean.mat
rm */02_preproc/*_clean_odd_full_ft.mat

# Dangerous - still part of pipeline?
#   SOLUTION: This was only for EEG12 because Sheila was in the process of re-cleaning.
#   These files will be overwritten once she's done anyways.
# ls -lh */03_events/*_behav_eeg_full_ft_02a.mat
#   SBJ02a saves SBJ_behav_eeg_full_ft_02a.mat
#   SBJ02c then loads that and converts to SBJ_behav_eeg_full_ft_final.mat
# ls -lh */03_events/*_behav_odd_full_ft_02a.mat
#   ODD02a saves SBJ_behav_odd_full_ft_02a.mat
#   ODD02c then loads that and converts to SBJ_behav_odd_full_ft_final.mat

rm */04_proc/*component*png
#   many component .png plots in EP SBJ and EEG04
rm -r */04_proc/ERP_stacks/
#   EP07 and EP09 have these directories with plots

# Random Files (I can't find them...)
## rm */02_preproc/*_eeg_full_ft_02a_tt.mat
## rm */02_preproc/*_odd_full_ft_02a_tt.mat
## rm EEG12_odd_full_ft_02a_odd.mat

## cd /data/EEG02/03_events/
## rm EEG02_behav_eeg_full_ft_clean.mat
## rm EEG02_behav_odd_full_ft_clean.mat
## rm EEG02_behav02a_eeg_full_ft_clean.mat
## rm EEG02_behav02b_odd_full_ft_clean.mat
## rm EEG02_behav02b_eeg_full_ft_clean.mat
## rm EEG02_behav02a_odd_full_ft_clean.mat
## cd ..
## cd /02_preproc/
## rm EEG02_clean_odd_full_ft.mat
## cd ..
## cd ..
## cd /EEG03/03_events/
## rm EEG03_behav_eeg_full_ft_clean.mat
## rm EEG03_behav_odd_full_ft_clean.mat
## rm EEG03_behav02a_eeg_full_ft_clean.mat
## rm EEG03_behav02b_odd_full_ft_clean.mat
## rm EEG02_behav02b_eeg_full_ft_clean.mat
## rm EEG03_behav02a_odd_full_ft_clean.mat
## cd ..
## cd /02_preproc/
## rm EEG03_clean_odd_full_ft.mat
## cd ..
## cd ..
## cd /EEG04/03_events/
## rm EEG04_behav_eeg_full_ft_clean.mat
## rm EEG04_behav_odd_full_ft_clean.mat
## rm EEG04_behav_oddball_eeg_full_ft_clean.mat
## rm EEG04_behav02a_eeg_full_ft_clean.mat
## rm EEG04_behav02b_eeg_full_ft_clean.mat
## rm EEG04_behav02b_odd_full_ft_clean.mat
## rm EEG04_behav02a_odd_full_ft_clean.mat
## cd ..
## cd ..
## cd /EEG05/02_preproc/
## rm EEG05_eeg_full_ft_02a_tt.mat
## rm EEG05_odd_full_ft_02a_odd.mat
## cd ..
## cd ..
## cd /EEG06/03_events/
## rm EEG06_behav_eeg_full_ft_clean.mat
## rm EEG06_behav_odd_full_ft_clean.mat
## rm EEG06_behav_oddball_eeg_full_ft_clean.mat
## rm EEG06_behav02a_eeg_full_ft_clean.mat
## rm EEG06_behav02b_eeg_full_ft_clean.mat
## rm EEG06_behav02b_odd_full_ft_clean.mat
## rm EEG06_behav02a_odd_full_ft_clean.mat
## cd ..
## cd ..
## cd /EEG07/03_events/
## rm EEG07_behav_eeg_full_ft_clean.mat
## rm EEG07_behav_odd_full_ft_clean.mat
## rm EEG07_behav02a_eeg_full_ft_clean.mat
## rm EEG07_behav02b_eeg_full_ft_clean.mat
## rm EEG07_behav02b_odd_full_ft_clean.mat
## rm EEG07_behav02a_odd_full_ft_clean.mat
## cd ..
## cd ..
## cd /EEG08/03_events/
## rm EEG08_behav_eeg_full_ft_clean.mat
## rm EEG07_behav_odd_full_ft_clean.mat
## rm EEG08_behav02a_eeg_full_ft_clean.mat
## rm EEG08_behav02b_eeg_full_ft_clean.mat
## rm EEG08_behav02b_odd_full_ft_clean.mat
## rm EEG08_behav02a_odd_full_ft_clean.mat
## cd ..
## cd ..
## cd /EEG10/03_events/
## rm EEG10_behav_eeg_full_ft_clean.mat
## rm EEG10_behav_odd_full_ft_clean.mat
## rm EEG10_behav02a_eeg_full_ft_clean.mat
## rm EEG10_behav02b_eeg_full_ft_clean.mat
## rm EEG10_behav02b_odd_full_ft_clean.mat
## rm EEG10_behav02a_odd_full_ft_clean.mat
## cd ..
## cd ..
## cd /EEG12/03_events/
## rm EEG12_behav_eeg_full_ft_02a.mat
## rm EEG12_behav_odd_full_ft_02a.mat
## rm EEG12_odd_full_ft_02a_odd.mat
## cd ..
## cd ..
## cd /EP06/03_events/
## rm EP06_behav_eeg_full_ft_clean.mat
## cd ..
## cd ..
## cd /EP07/04_preproc/
## rm EP07component*
## rm -r ERP_stacks/
## cd ..
## cd ..
## cd /EP09/04_preproc/
## rm -r ERP_Stacks/
## rm EP09component*
## cd ..
## cd ..
## cd /EP08/04_preproc/
## rm EP08component*
## cd ..
## cd ..
## cd /EP10/04_preproc/
## rm EP10component*
## cd ..
## cd ..
## cd /EP11/04_preproc/
## rm EP11component*
 


