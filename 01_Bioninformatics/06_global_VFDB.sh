#!/usr/bin/env bash
## Last update 04-09-2025 by Suleiman & AbdulAziz
#SBATCH --job-name=global_vfdb
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=36:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/slurm-%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/slurm-%j.err


# ---- Activate env and PATH ----
eval "$(conda shell.bash hook)" || true
export PATH="$HOME/anaconda3/bin:$PATH"

# ---- Paths ----
BASE_PATH="/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/DATA"
VFDB_DB="/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Shotgun-metagenomics/VFDB/setA"
DIAMOND_path="/srv/lustre01/project/mmrd-cp3fk69sfrq/shared/diamond"
GLOBAL_DEREP="${BASE_PATH}/Cancer_metagenomics_global_derep"
REP_GENOMES="${GLOBAL_DEREP}/output/dereplicated_genomes"
VFDB_OUT="${GLOBAL_DEREP}/VFDB_OUT"
THREADS=${SLURM_CPUS_PER_TASK:-56}
MEM_GB=96

mkdir -p "$VFDB_OUT/annotations" "$VFDB_OUT/hit_fastas" "$VFDB_OUT/indices" "$VFDB_OUT/mappings" "$VFDB_OUT/abundance"

# ---- Runtime log setup ----
RUNTIME_LOG="$GLOBAL_DEREP/runtime_metrics_global_VFDB3.csv"
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
# Step 3: Build Bowtie2 index
############################################
#echo "[Step 3] Building Bowtie2 index"
#conda activate bowtie2 || exit 1
#step_start=$(date +%s)
#bowtie2-build "$VFDB_OUT/combined_vfdb_ref.fna" "$VFDB_OUT/indices/VFDB_index"
#step_end=$(date +%s)
#log_step "global" "Bowtie2_index" "$step_start" "$step_end" "$VFDB_OUT/combined_vfdb_ref.fna" "$VFDB_OUT/indices/VFDB_index"
#echo "  Bowtie2 index built"


############################################
# Step 4: Map reads to VFDB
############################################
eecho "[Step 4] Mapping reads"
conda activate bowtie2 || exit 1

for PROJECT in PRJEB72526; do
    SRA_PATH="$BASE_PATH/$PROJECT"

    for sample in $(cat "$SRA_PATH/SRA_ids3"); do
        R1="$SRA_PATH/$sample/${sample}_R1.decontam.paired.fq.gz"
        R2="$SRA_PATH/$sample/${sample}_R2.decontam.paired.fq.gz"
        SAM_OUT="$VFDB_OUT/mappings/${sample}_VFDB.sam"

        # Skip if SAM file already exists
        if [ -f "$SAM_OUT" ]; then
            echo " ⏩ Skipping $sample (SAM already exists)"
            continue
        fi

        step_start=$(date +%s)

        bowtie2 -x "$VFDB_OUT/indices/VFDB_index" \
                -1 "$R1" \
                -2 "$R2" \
                -S "$SAM_OUT" \
                -p "$THREADS" \
                --quiet

        step_end=$(date +%s)

        log_step "$sample" "Bowtie2_map" "$step_start" "$step_end" "$R1" "$SAM_OUT"
        echo " ✔ Mapping done for $sample"
    done
done


############################################
# Step 4.5: Convert SAM to sorted BAM
############################################
#echo "[Step 4.5] Converting SAM to BAM"
#conda activate samtools || exit 1
#for sam in "$VFDB_OUT/mappings"/*.sam; do
#sample=$(basename "$sam" _VFDB.sam)
#BAM="$VFDB_OUT/mappings/${sample}_VFDB.sorted.bam"
#step_start=$(date +%s)
#samtools view -bS "$sam" | samtools sort -@ "$THREADS" -o "$BAM"
#samtools index "$BAM"
#step_end=$(date +%s)
#log_step "$sample" "SAM_to_BAM" "$step_start" "$step_end" "$sam" "$BAM"
#echo " ✔ Sorted BAM complete for $sample"
#done


############################################
# Step 5: CoverM abundance
############################################
#echo "[Step 5] Quantifying VFDB abundance"
#conda activate coverm || exit 1
#for bam in "$VFDB_OUT/mappings"/*_VFDB.sorted.bam; do
#sample=$(basename "$bam" _VFDB.sorted.bam)
#step_start=$(date +%s)
#coverm contig \
#--bam-files "$bam" \
#--methods rpkm tpm count \
#--min-read-percent-identity 80 \
#--min-read-aligned-percent 80 \
#--threads "$THREADS" \
#--output-file "$VFDB_OUT/abundance/${sample}_VFDB_abundance.tsv"
#step_end=$(date +%s)
#log_step "$sample" "CoverM_VFDB_abundance" "$step_start" "$step_end" "$bam" "$VFDB_OUT/abundance/${sample}_VFDB_abundance.tsv"
#echo " ✔ CoverM complete for $sample"
#done


#echo "🏁 Pipeline complete: Global VFDB abundance quantification ready."
