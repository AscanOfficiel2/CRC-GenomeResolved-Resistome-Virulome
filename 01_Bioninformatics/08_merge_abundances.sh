#!/usr/bin/env bash
## Script 2: Merge cleaned files + extract GTDB genome-species maps
## Last update: FINAL VERSION
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=35:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/merge-%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/merge-%j.err


# ---- Base Configuration ----
BASE_DIR="/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/DATA/Cancer_metagenomics_global_derep"
MERGED_DIR="$BASE_DIR/merged_outputs"
mkdir -p "$MERGED_DIR"

# Cleaned abundance paths
declare -A INPUT_PATHS=(
  [gtdb]="$BASE_DIR/coverm_cleaned"
  [arg]="$BASE_DIR/arg_abundance/abundance_cleaned"
  [vfdb]="$BASE_DIR/VFDB_OUT/abundance_cleaned"
)

# ---- Step 1: Merge cleaned abundance files ----
for TYPE in gtdb arg vfdb; do
  echo "[Merging cleaned $TYPE abundance files...]"
  CLEAN_DIR="${INPUT_PATHS[$TYPE]}"
  OUTFILE="$MERGED_DIR/${TYPE}_combined.tsv"
  HEADER_WRITTEN=false

  for f in "$CLEAN_DIR"/*.tsv; do
    [[ -e "$f" ]] || continue

    if [[ "$HEADER_WRITTEN" == false ]]; then
      head -n1 "$f" > "$OUTFILE"
      HEADER_WRITTEN=true
    fi

    tail -n +2 "$f" >> "$OUTFILE"
  done

  echo "  ✔ Combined file written: $OUTFILE"
done

# ---- Step 2: Extract genome-to-species maps ----
echo "[Extracting GTDB genome-to-species mappings...]"

GTDB_SUMMARY_BAC="$BASE_DIR/taxonomy/gtdbtk/classify/gtdbtk.bac120.summary.tsv"
GTDB_SUMMARY_ARCH="$BASE_DIR/taxonomy/gtdbtk/classify/gtdbtk.ar53.summary.tsv"

MAP_BAC="$MERGED_DIR/gtdb_genome_species.bac120_map.tsv"
MAP_ARCH="$MERGED_DIR/gtdb_genome_species.ar53_map.tsv"

awk -F '\t' 'NR > 1 {print $1 "\t" $2}' "$GTDB_SUMMARY_BAC" > "$MAP_BAC"
awk -F '\t' 'NR > 1 {print $1 "\t" $2}' "$GTDB_SUMMARY_ARCH" > "$MAP_ARCH"

echo "  ✔ Genome-to-species maps saved:"
echo "     - $MAP_BAC"
echo "     - $MAP_ARCH"


echo "✅ All merged tables and species maps are ready in: $MERGED_DIR"
