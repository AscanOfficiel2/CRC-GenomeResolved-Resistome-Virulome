#!/bin/bash

MAIN_DIR="/home/rachidlab/Downloads/MAGS/projects/MAGS_QUALITY"
LIST_FILE="$MAIN_DIR/list"

cd "$MAIN_DIR" || exit

# --- Setup Directories ---
mkdir -p collected_checkm2
mkdir -p collected_fastqscreen_before1
mkdir -p collected_fastqscreen_before2
# New directories for Quast files
mkdir -p collected_filtered_quast
mkdir -p collected_simplified_quast

echo "Reading folders from list..."
while IFS= read -r folder; do
    echo "Processing: $folder"
    # Full path to the folder
    FOLDER_PATH="$MAIN_DIR/$folder"

    if [ -d "$FOLDER_PATH" ]; then
        # Collect CheckM2 files
        find "$FOLDER_PATH" -type f -name "*_checkm2.csv" -exec cp {} collected_checkm2/ \;

        # Collect Fastqscreen_before1 files
        find "$FOLDER_PATH" -type f -name "*_Fastqscreen_before1.csv" -exec cp {} collected_fastqscreen_before1/ \;

        # Collect Fastqscreen_before2 files
        find "$FOLDER_PATH" -type f -name "*_Fastqscreen_before2.csv" -exec cp {} collected_fastqscreen_before2/ \;

        # Collect Filtered Quast files (new)
        find "$FOLDER_PATH" -type f -name "*_filtered_quast.csv" -exec cp {} collected_filtered_quast/ \;

        # Collect Simplified Quast files (new)
        find "$FOLDER_PATH" -type f -name "*_simplified_quast.csv" -exec cp {} collected_simplified_quast/ \;
    else
        echo "Warning: Folder not found → $folder"
    fi
done < "$LIST_FILE"

# --- Function to concatenate CSV files with header management ---
concatenate_csvs() {
    local input_dir=$1
    local output_file=$2
    local first=1

    echo "Concatenating all files from $input_dir..."
    # Clear or create the output file
    > "$output_file"

    # Check if any files exist in the directory before proceeding
    shopt -s nullglob
    files=("$input_dir"/*.csv)
    shopt -u nullglob

    if [ ${#files[@]} -eq 0 ]; then
        echo "No files found in $input_dir. Skipping concatenation for this set."
        return
    fi

    for file in "${files[@]}"; do
        if [ $first -eq 1 ]; then
            cat "$file" >> "$output_file"
            first=0
        else
            # Append content starting from the second line (excluding header)
            tail -n +2 "$file" >> "$output_file"
        fi
    done
    echo "Done!"
    echo "Combined file created: $MAIN_DIR/$output_file"
}

# --- Execute Concatenation for all five sets ---
echo "--- Processing CheckM2 files ---"
concatenate_csvs collected_checkm2 all_checkm2_combined.csv
echo ""

echo "--- Processing Fastqscreen Before 1 files ---"
concatenate_csvs collected_fastqscreen_before1 all_Fastqscreen_before1_combined.csv
echo ""

echo "--- Processing Fastqscreen Before 2 files ---"
concatenate_csvs collected_fastqscreen_before2 all_Fastqscreen_before2_combined.csv
echo ""

# Execute Concatenation for Quast sets (new)
echo "--- Processing Filtered Quast files ---"
concatenate_csvs collected_filtered_quast all_filtered_quast_combined.csv
echo ""

echo "--- Processing Simplified Quast files ---"
concatenate_csvs collected_simplified_quast all_simplified_quast_combined.csv
echo ""

echo "All tasks complete!"

