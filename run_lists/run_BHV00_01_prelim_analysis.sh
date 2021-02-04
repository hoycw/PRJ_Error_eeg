#!/bin/sh

# Move to scripts directory 
root_dir="/Volumes/hoycw_clust/PRJ_Error_eeg/"
cd $root_dir/scripts/

# Load the correct SBJ list
echo "Please enter the desired SBJ list:"
read list_name

declare -a SBJs
SBJs=(`cat "${root_dir}scripts/SBJ_lists/${list_name}.sbj"`)
#    =("P2" "P3" "P7" "P1" "P6" "P4" "P5" "P8" "P9"\
#        "IR75" "IR76" "colin_vec" "IR79" "colin_novec" "IR84" "CP23" "CP22" "IR82" "Rana_1.6"\
#        "IR71" "CP241" "Adi_1.7" "IR63" "IR62" "IR72" "IR60" "IR67" "IR66" "IR65" "IR77" "IR78"\
#        "IR57" "IR69" "IR68" "CP242" "IR74")

# Run BHV scripts for these SBJs
for sbj in "${SBJs[@]}"; do
    echo "=================================================================\n"
    echo "Running ${sbj}\n"
    echo "=================================================================\n"
    python BHV00_extract.py ${sbj}
#    if [ ! -f ${root_dir}/data/${sbj}/03_events/${sbj}_prdm_vars.pkl ]; then
#        echo "pkl of paradigm vars not found, running for ${sbj}..."
#        python BHV00_extract.py ${sbj}
#    else
#        echo "pkl found, continuing to BHV01 for ${sbj}..."
#    fi
#    python BHV01_prelim_analysis.py ${sbj}
done

