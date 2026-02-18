#!/usr/bin/env bash
## Script 3: Annotate GTDB, ARG and VFDB with species info
## EXACT ORIGINAL LOGIC – ONLY PATHS FIXED + MINIMAL GENOME TRIM FIX
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=35:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/slurm-%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/slurm-%j.err

# =========================================
#        CONFIGURATION
# =========================================
BASE_DIR="/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/DATA/Cancer_metagenomics_global_derep"
MERGED_DIR="$BASE_DIR/merged_outputs"

GTDB_IN="$MERGED_DIR/gtdb_combined.tsv"
ARG_IN="$MERGED_DIR/arg_combined.tsv"
VFDB_IN="$MERGED_DIR/vfdb_combined.tsv"

GTDB_MAP_BAC="$MERGED_DIR/gtdb_genome_species.bac120_map.tsv"
GTDB_MAP_ARCH="$MERGED_DIR/gtdb_genome_species.ar53_map.tsv"

FINAL_OUT="$MERGED_DIR/final_outputs"
mkdir -p "$FINAL_OUT"

# =========================================
# STEP 1: Build unified genome-to-species map
# =========================================
echo "[1] Building genome-to-species map..."

awk 'FNR==1 && NR!=1{next}{print}' \
    "$GTDB_MAP_BAC" "$GTDB_MAP_ARCH" \
    > "$FINAL_OUT/genome_to_species.tsv"

echo "   → genome_to_species.tsv created."


# =========================================
# STEP 2: Annotate GTDB (Genome is column 3 — EXACT ORIGINAL LOGIC)
# =========================================
echo "[2] Annotating GTDB table with species..."

awk -v map="$FINAL_OUT/genome_to_species.tsv" '
BEGIN {
    FS=OFS="\t"
    while ((getline < map) > 0) {
        species[$1] = $2
    }
}
NR==1 {
    print $0, "Species"
    next
}
{
    genome = $3
    if (genome in species)
         print $0, species[genome]
    else print $0, "Unmapped"
}
' "$GTDB_IN" > "$FINAL_OUT/gtdb_with_species.tsv"

echo "   → gtdb_with_species.tsv created."


# =========================================
# STEP 3: Extract genome IDs from ARG and VFDB (EXACT ORIGINAL LOGIC)
# =========================================

extract_arg() {
    infile=$1
    outfile=$2

    echo "[3A] Extracting genome IDs from ARG..."

    awk 'BEGIN{FS=OFS="\t"}
    NR==1 { print $0,"Genome"; next }
    {
        header=$2
        genome="Unmapped"
        if (match(header,/sample=([^|]+)/,m))
            genome=m[1]
        print $0, genome
    }' "$infile" > "$outfile"
}

extract_vfdb() {
    infile=$1
    outfile=$2

    echo "[3B] Extracting genome IDs from VFDB..."

    awk 'BEGIN{FS=OFS="\t"}
    NR==1 { print $0,"Genome"; next }
    {
        header=$2
        genome=header
        sub(/\|.*/, "", genome)
        print $0, genome
    }' "$infile" > "$outfile"
}

extract_arg   "$ARG_IN"  "$FINAL_OUT/arg_with_genome.tsv"
extract_vfdb  "$VFDB_IN" "$FINAL_OUT/vfdb_with_genome.tsv"

echo "   → genome extraction complete."


# =========================================
# STEP 4: Annotate ARG + VFDB with species
# FIX APPLIED: trim spaces around genome IDs
# =========================================

annotate_table() {
    infile=$1
    outfile=$2

    awk -v map="$FINAL_OUT/genome_to_species.tsv" '
    BEGIN {
        FS=OFS="\t"
        while ((getline line < map) > 0) {
            gsub(/^[ \t]+|[ \t]+$/, "", line)
            split(line, a, "\t")
            gsub(/^[ \t]+|[ \t]+$/, "", a[1])
            species[a[1]] = a[2]
        }
    }

    NR==1 {
        print $0, "Species"
        next
    }

    {
        genome=$(NF)
        gsub(/^[ \t]+|[ \t]+$/, "", genome)

        if (genome in species)
             print $0, species[genome]
        else print $0, "Unmapped"
    }
    ' "$infile" > "$outfile"
}

echo "[4] Annotating ARG + VFDB with species..."

annotate_table "$FINAL_OUT/arg_with_genome.tsv"  "$FINAL_OUT/arg_with_species.tsv"
annotate_table "$FINAL_OUT/vfdb_with_genome.tsv" "$FINAL_OUT/vfdb_with_species.tsv"

echo "   → arg_with_species.tsv"
echo "   → vfdb_with_species.tsv"


# =========================================
# DONE
# =========================================
echo "=============================================="
echo "   ALL DONE — Outputs in: $FINAL_OUT"
echo "=============================================="

# =========================================
# STEP 2B: Separate mapping for ARCHAEA ONLY
# =========================================
echo "[2B] Performing separate annotation for archaeal genomes..."

ARCH_OUT="$FINAL_OUT/gtdb_arch_with_species.tsv"

awk -v map="$GTDB_MAP_ARCH" '
BEGIN {
    FS=OFS="\t"
    # Load archaeal genome→species map
    while ((getline < map) > 0) {
        arch_species[$1] = $2
    }
}
NR==1 {
    header = $0
    print header, "Species"
    next
}
{
    genome = $3

    # Check only archaeal genomes
    if (genome in arch_species) {
        print $0, arch_species[genome]
    }
}
' "$GTDB_IN" > "$ARCH_OUT"

echo "   → archaeal species annotation saved to: $ARCH_OUT"
