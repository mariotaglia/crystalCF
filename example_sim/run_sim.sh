#!/bin/bash

# List to check.
file_list=("F_tot_gcanon.dat" "F_trans.dat" "F_trans_sv.dat" "F_vdW.dat" "F_conf.dat" "F_conf_sv.dat" "F_mixs.dat" "F_HS.dat")

output_file="missing_paths.txt"
echo "" > "$output_file"

# Search all dir without subfolders with a DEFINITIONS file.
find . -type d | while read -r dir; do
    # Check if DEFINTIONS.txt exists.
    if [[ -f "$dir/DEFINITIONS.txt" ]]; then
        if [[ -z $(find "$dir" -mindepth 1 -type d) ]]; then
            found_issue=false
            for filename in "${file_list[@]}"; do
                file="$dir/$filename"
                # Checks if the file is empty or doesn't exist.
                if [[ ! -f "$file" || ! -s "$file" ]]; then
                    if [[ "$found_issue" == false ]]; then
                        echo "$dir" >> "$output_file"
                        found_issue=true
                    fi
                fi
            done
        fi
    fi
done

# Execute the sbatch command to rerun each missing path.
sort -u "$output_file" | while read -r folder; do
    if [[ "$folder" != "$(pwd)" ]]; then
        (cd "$folder" && sbatch tosubmit.sh)
    fi
done
