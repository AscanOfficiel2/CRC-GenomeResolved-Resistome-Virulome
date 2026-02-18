#!/usr/bin/env bash
## Script 1: Clean abundance files (GTDB, ARG, VFDB)
## Last update: FINAL VERSION
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=35:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/clean-%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/clean-%j.err

# ---- Base Configuration ----
BASE_DIR="/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/DATA/Cancer_metagenomics_global_derep"

# Tool input folders
TOOLS=("coverm" "arg_abundance/abundance" "VFDB_OUT/abundance")
TOOL_LABELS=("GTDB" "ARG" "VFDB")

# ---- Processing Loop ----
for i in "${!TOOLS[@]}"; do
  TOOL_PATH="${TOOLS[$i]}"
  TOOL_LABEL="${TOOL_LABELS[$i]}"

  INPUT_DIR="$BASE_DIR/$TOOL_PATH"
  OUTPUT_DIR="$BASE_DIR/${TOOL_PATH}_cleaned"
  mkdir -p "$OUTPUT_DIR"

  echo "[Processing $TOOL_LABEL] Input: $INPUT_DIR → Output: $OUTPUT_DIR"

  for f in "$INPUT_DIR"/*.tsv; do
    [[ -e "$f" ]] || continue
    sample=$(basename "$f" .tsv | sed 's/_\(arg\|VFDB\)_abundance//')
    outfile="$OUTPUT_DIR/${sample}_${TOOL_LABEL}_cleaned.tsv"

    echo "  - Cleaning sample: $sample"

    # Copy header
    head -n 1 "$f" > "$outfile"

    # Filter non-zero rows and append with sample ID
    awk -v sample="$sample" 'NR>1 {
      keep = 0
      for (i = 2; i <= NF; i++) {
        if ($i ~ /^[0-9.]+$/ && $i+0 != 0) {
          keep = 1
          break
        }
      }
      if (keep) {
        print sample "\t" $0
      }
    }' "$f" >> "$outfile"
  done
done

echo "✅ All abundance tables cleaned and saved to *_cleaned folders."
