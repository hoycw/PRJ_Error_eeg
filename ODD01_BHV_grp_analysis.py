if [ -d "/Volumes/hoycw_clust/PRJ_Error_eeg/" ]; then
    root_dir="/Volumes/hoycw_clust/PRJ_Error_eeg/"
else
    root_dir="/Users/sheilasteiner/Desktop/Knight_Lab/PRJ_Error_eeg/"
fi
cd $root_dir/scripts/

# Load the correct SBJ list
echo "Please enter the desired SBJ list:"
read list_name

declare -a SBJs
SBJs=(`cat "${root_dir}scripts/SBJ_lists/${list_name}.sbj"`)
for sbj in "${SBJs[@]}"; do
    echo "=================================================================\n"
    echo "Running ${sbj}\n"
    echo "=================================================================\n"
    
    python ODD01_grp_analysis.py ${sbj}
done

