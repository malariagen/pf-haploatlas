#!/bin/bash

# This script will be run after generate_gene_summary.py to curate the logs. 

# Today's date in the format YYYY-MM-DD
today_date=$(date +"%Y-%m-%d")

# Archive and remove the old job_logs.json
destination_job_logs="job_logs_archive"
old_job_log="job_logs.json"

# Check if old job_logs.json exists
if [ ! -e "$old_job_log" ]; then
    echo "Error: $old_job_log not found."
    exit 1
fi

# Move old job_logs.json to the destination directory
mv "$old_job_log" "$destination_job_logs/${today_date}_${old_job_log}"

# Check if mv command executed successfully
if [ $? -ne 0 ]; then
    echo "Error: Failed to move $old_job_log to $destination_job_logs."
    exit 1
else
    echo "Successfully moved $old_job_log to $destination_job_logs."
fi

#### PARSE THE OUTPUTS

# Define a variable to store the sum of runtimes
total_runtime_seconds=0

# Loop through all output files starting with "output."
for file in output.*; do
    # Check if the file exists and is not empty
    if [ -s "$file" ]; then
        
        # Extract job ID
        job_id=$(grep -oP 'Job \K\d+' "$file")
        # Extract runtime using grep and awk
        runtime_seconds=$(grep -oP 'Run time : +\K\d+' "$file")
        
        # Add runtime to the total
        total_runtime_seconds=$((total_runtime_seconds + runtime_seconds))
    fi
done

# Convert total runtime to hours
total_runtime_hours=$(awk "BEGIN {printf \"%.2f\", $total_runtime_seconds / 216000}")

# Print the job ID 
echo "Job ID: $job_id"

# Print the total runtime in hours
echo "Total Runtime: $total_runtime_hours hours"

# Remove files starting with "error." and "output."
rm -f error.* output.*

# Update the JSON file with job_id and total_runtime_hours
jq --arg job_id "$job_id" --arg total_runtime_hours "$total_runtime_hours" \
   '.job_id = $job_id | .job_runtime = $total_runtime_hours' \
   gene_logs.json > job_logs.json
   
echo "Logs stored in job_logs.json"

#### ARCHIVE THE SCRIPT 

# Your existing script file
script_file="generate_gene_summary.py"

# Destination folder
destination_folder="script_archive"

# Destination path with the new filename
destination_path="$destination_folder/${today_date}_${script_file%.py}.py"

# Copy the file to the destination folder and rename it
cp "$script_file" "$destination_path"

echo "Script copied to: $destination_path"

#### MOVE PKL FILES

# # Source folder where your pkl.xz files are located
# source_folder="."

# # Destination folder for pkl files
# destination_folder="pkl_files"

# # Create the destination folder if it doesn't exist
# mkdir -p "$destination_folder"

# # Move all files ending with "pkl.xz" to the destination folder
# mv "$source_folder"/*.pkl.xz "$destination_folder/"

# echo "Files moved to: $destination_folder"

# Source folder where your pkl.xz files are located
source_folder="out_pkls"

# Destination folder for pkl files
destination_folder=../../app/files/"${today_date}_pkl_files"

# Create the destination folder if it doesn't exist
mkdir -p "$destination_folder"

# Move all files ending with "pkl.xz" to the destination folder
# First, copy the files
scp "$source_folder"/*.pkl.xz "$destination_folder/"

# Check if the copy was successful
if [ $? -eq 0 ]; then
    # Count the number of files in the destination folder
    file_count=$(ls -1 "$destination_folder"/*.pkl.xz 2>/dev/null | wc -l)

    # Check if the file count is 5102
    if [ "$file_count" -eq 5102 ]; then
        # If there are exactly 5102 files, delete the files from the source folder
        rm "$source_folder"/*.pkl.xz
        rmdir ${source_fol}
        echo "Files moved to: $destination_folder and deleted from source."
    else
        echo "File copy succeeded, but the destination folder does not contain exactly 5102 files. No files were deleted from the source."
    fi
else
    echo "File copy failed. No files were deleted from the source."
fi