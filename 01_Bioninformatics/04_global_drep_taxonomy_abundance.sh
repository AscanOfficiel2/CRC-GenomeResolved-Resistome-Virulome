#!/usr/bin/env bash
## Global ARG Abundance Quantification from Dereplicated MAGs
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=35:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/slurm-%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/slurm-%j.err

# -------------------------------
# Activate conda and set env
# -------------------------------
eval "$(conda shell.bash hook)" || true
export PATH="$HOME/anaconda3/bin:$PATH"

# -------------------------------
# Paths and settings
# -------------------------------
BASE_PATH="/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/DATA"
CHECKM2DB="/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Shotgun-metagenomics/checkm2_databases/CheckM2_database/"
GTDBTK_DB="/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Shotgun-metagenomics/gtbk_db_new/release226/"
THREADS=${SLURM_CPUS_PER_TASK:-56}
MEM_GB=96

# -------------------------------
# Output dirs
# -------------------------------
GLOBAL_OUT="$BASE_PATH/Cancer_metagenomics_global_derep"
ALL_GENOMES="$GLOBAL_OUT/genomes"
GTDB_OUT="$GLOBAL_OUT/taxonomy/gtdbtk"
COVERM_OUT="$GLOBAL_OUT/coverm"
REP_GENOMES="$GLOBAL_OUT/output/dereplicated_genomes"
mkdir -p "$ALL_GENOMES" "$GTDB_OUT" "$COVERM_OUT" "$REP_GENOMES"

# -------------------------------
# Runtime log setup
# -------------------------------
RUNTIME_LOG="$GLOBAL_OUT/runtime_metrics_global_new1.csv"
if [[ ! -f "$RUNTIME_LOG" ]]; then
  echo "Sample,Step,StartTime,EndTime,Duration_sec,Threads,Mem_GB,InputSize,OutputSize" > "$RUNTIME_LOG"
fi

log_step () {
  local sample=$1 step=$2 start=$3 end=$4 infile=$5 outfile=$6
  local dur=$((end - start))
  local insize=$( [ -f "$infile" ] && du -sh "$infile" | cut -f1 || echo "NA" )
  local outsize=$( [ -f "$outfile" ] && du -sh "$outfile" | cut -f1 || echo "NA" )
  echo "$sample,$step,$start,$end,$dur,$THREADS,$MEM_GB,$insize,$outsize" >> "$RUNTIME_LOG"
}

############################################
# Step 1: Collect HQ MAGs across projects
############################################
#echo "[Step 1] Collecting HQ MAGs for global dereplication"
#for PROJECT in PRJDB4176 PRJEB6070 PRJEB7774 PRJEB27928 PRJEB72525 PRJEB72526 PRJNA389927 PRJNA1167935; do
 # SRA_PATH="$BASE_PATH/$PROJECT"
  #if [[ ! -d "$SRA_PATH" ]]; then
   # echo "⚠️  Project path not found: $SRA_PATH"
    #continue
  #fi

  #while read -r SAMPLE; do
   # [[ -z "$SAMPLE" ]] && continue
    #CHECKM_TSV="$SRA_PATH/$SAMPLE/checkm2_$SAMPLE/quality_report.tsv"
   # BIN_DIR="$SRA_PATH/$SAMPLE/bins"

    #if [[ -f "$CHECKM_TSV" ]]; then
     # echo "  ➤ Checking $SAMPLE in $PROJECT"
      #awk -F'\t' 'NR>1 && $2 >= 80 && $3 <= 10 {print $1}' "$CHECKM_TSV" | while read -r bin; do
       # src="$BIN_DIR/${bin}.fa"
        #dest="$ALL_GENOMES/${PROJECT}_${SAMPLE}_${bin}.fa"
        #if [[ -f "$src" ]]; then
         # cp "$src" "$dest"
          #echo "    Copied: $src → $dest"
        #fi
      #done
   # fi
  #done < "$SRA_PATH/SRA_ids"
#done

############################################
# Step 2: dRep dereplication (safe input list)
############################################
#echo "[Step 2] Running dRep dereplication"
#conda activate drep || exit 1

# Build a file list instead of wildcards
#find "$ALL_GENOMES" -type f -name "*.fa" > "$GLOBAL_OUT/genome_list.txt"
#GENOME_COUNT=$(wc -l < "$GLOBAL_OUT/genome_list.txt")
#echo "  ➤ Found $GENOME_COUNT genome files for dereplication"

#step_start=$(date +%s)
#dRep dereplicate "$GLOBAL_OUT/output" \
 # -g "$GLOBAL_OUT/genome_list.txt" \
  #-p "$THREADS" \
  #--ignoreGenomeQuality \
  #-sa 0.95 -nc 0.30 -pa 0.90 --S_algorithm fastANI
#step_end=$(date +%s)

#log_step "global" "dRep" "$step_start" "$step_end" "$ALL_GENOMES" "$GLOBAL_OUT/output"
#echo "✅ dRep finished in $((step_end - step_start)) seconds"

# Convert .fa → .fna (for CoverM)
#if compgen -G "$GLOBAL_OUT/output/dereplicated_genomes/*.fa" > /dev/null; then
  #for f in "$GLOBAL_OUT/output/dereplicated_genomes/"*.fa; do
   # cp "$f" "${f%.fa}.fna"
  #done
  #echo "✅ Converted dereplicated genomes from .fa → .fna"
#else
 # echo "⚠️  No dereplicated genomes found — skipping .fa → .fna conversion."
#fi

############################################
# Step 3: GTDB-Tk taxonomy assignment
############################################
#echo "[Step 3] Running GTDB-Tk"
#conda activate gtdbtk-2.4.1 || exit 1
#export GTDBTK_DATA_PATH="$GTDBTK_DB"

#step_start=$(date +%s)
#gtdbtk classify_wf \
 # --genome_dir "$REP_GENOMES" \
  #--out_dir "$GTDB_OUT" \
  #--cpus "$THREADS" \
  #--extension fa \
  #--skip_ani_screen
#step_end=$(date +%s)

#log_step "global" "GTDB-Tk" "$step_start" "$step_end" "$REP_GENOMES" "$GTDB_OUT"
#echo "✅ GTDB-Tk finished in $((step_end - step_start)) seconds"

############################################
# Step 4: CoverM abundance profiling (optional)
############################################
echo "[Step 4] Running CoverM on all samples"
conda activate coverm || exit 1
for PROJECT in PRJDB4176; do
  SRA_PATH="$BASE_PATH/$PROJECT"
  while read -r SAMPLE; do
    [[ -z "$SAMPLE" ]] && continue
    SAMPLE_DIR="$SRA_PATH/$SAMPLE"
    step_start=$(date +%s)
    coverm genome \
      -1 "$SAMPLE_DIR/${SAMPLE}_R1.decontam.paired.fq.gz" \
      -2 "$SAMPLE_DIR/${SAMPLE}_R2.decontam.paired.fq.gz" \
      --genome-fasta-directory "$REP_GENOMES" \
      --threads "$THREADS" \
      --output-format sparse \
      --methods rpkm tpm covered_fraction relative_abundance trimmed_mean mean variance \
      --output-file "$COVERM_OUT/${SAMPLE}_coverm.tsv"
    step_end=$(date +%s)
    log_step "$SAMPLE" "CoverM" "$step_start" "$step_end" "$REP_GENOMES" "$COVERM_OUT/${SAMPLE}_coverm.tsv"
  done < "$SRA_PATH/SRA_ids"
done

echo "🎯 Global dereplication, taxonomy, and abundance workflow complete."
