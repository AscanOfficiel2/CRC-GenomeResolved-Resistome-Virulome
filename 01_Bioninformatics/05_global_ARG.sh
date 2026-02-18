#!/usr/bin/env bash
## Global ARG Abundance Quantification from Dereplicated MAGs
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=35:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/slurm-%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/slurm-%j.err

# ---- Environment ----
eval "$(conda shell.bash hook)" || true
export PATH="$HOME/anaconda3/bin:$PATH"


# ---- Input paths ----
BASE_PATH="/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/DATA"
PROJECT="Cancer_metagenomics_global_derep"
DREP_OUT="$BASE_PATH/Cancer_metagenomics_global_derep"
REP_GENOMES="$DREP_OUT/output/dereplicated_genomes"
ARG_DIR="$DREP_OUT/arg_abundance"
THREADS=${SLURM_CPUS_PER_TASK:-56}
MEM_GB=96

mkdir -p "$ARG_DIR/annotations" "$ARG_DIR/indices" "$ARG_DIR/mappings" "$ARG_DIR/abundance"

# ---- Runtime log ----
RUNTIME_LOG="$DREP_OUT/runtime_metrics_${PROJECT}_ARG2_sam_bam.csv"
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
# Step 1: Run RGI on dereplicated MAGs
############################################
#echo "[Step 1] Running RGI on MAG nucleotide contigs"
#conda activate rgi || exit 1
#for MAG in "$REP_GENOMES"/*.fna; do
#  base=$(basename "$MAG" .fna)
#  step_start=$(date +%s)
#  rgi main \
#    --input_sequence "$MAG" \
#    --output_file "$ARG_DIR/annotations/${base}_rgi" \
#    --input_type contig \
#    --num_threads "$THREADS" \
#    --clean \
#    --data wgs
#  step_end=$(date +%s)
#  log_step "$base" "RGI_annotation" "$step_start" "$step_end" "$MAG" "$ARG_DIR/annotations/${base}_rgi.txt"
#  echo "  RGI complete for $base"
#done

############################################
# Step 2: Extract ARG DNA sequences to FASTA
############################################
#echo "[Step 2] Extracting ARG DNA sequences from RGI TXT outputs"
#ARG_FASTA="$ARG_DIR/ARG_reference.fna"
#> "$ARG_FASTA"

#for txt in "$ARG_DIR/annotations"/*.txt; do
#  base=$(basename "$txt" _rgi.txt)
#  step_start=$(date +%s)
#  awk -F'\t' -v BASE="$base" 'NR > 1 && $18 != "" {
#    ORF_ID = $1; CONTIG = $2; START = $3; STOP = $4; ORIENT = $5; DNA = $18;
#    split($9, g, " "); GENE = g[1]
#    gsub(/[^a-zA-Z0-9_().-]/, "_", ORF_ID); gsub(/[^a-zA-Z0-9_().-]/, "_", CONTIG); gsub(/[^a-zA-Z0-9_().-]/, "_", GENE); gsub(/[^a-zA-Z0-9_().-]/, "_", BASE);
#    printf(">%s|%s|%s-%s|%s|gene=%s|sample=%s\n%s\n", CONTIG, ORF_ID, START, STOP, ORIENT, GENE, BASE, DNA)
#  }' "$txt" >> "$ARG_FASTA"
#  step_end=$(date +%s)
#  log_step "$base" "Extract_ARGs" "$step_start" "$step_end" "$txt" "$ARG_FASTA"
#done

#if [[ ! -s "$ARG_FASTA" ]]; then
#  echo "❌ ERROR: No ARGs extracted. Check RGI TXT outputs."
#  exit 1
#fi

#echo "✅ ARG reference FASTA created: $ARG_FASTA"

############################################
# Step 3: Build Bowtie2 index
############################################
#echo "[Step 3] Building Bowtie2 index"
#conda activate bowtie2 || exit 1
#step_start=$(date +%s)
#bowtie2-build "$ARG_FASTA" "$ARG_DIR/indices/arg_index"
#step_end=$(date +%s)
#log_step "$PROJECT" "Bowtie2_index" "$step_start" "$step_end" "$ARG_FASTA" "$ARG_DIR/indices/arg_index"
#echo "✅ Bowtie2 index built"



############################################
# Step 4: Map reads and generate BAM # run this first before sam and bam for all projects
############################################
# Step 4: Map reads and generate BAM
############################################
#echo "[Step 4] Mapping reads from all samples"

#for PROJECT_DIR in PRJEB6070; do
#  SRA_PATH="$BASE_PATH/$PROJECT_DIR"

#  for sample in $(cat "$SRA_PATH/SRA_ids1"); do
#    R1="$SRA_PATH/$sample/${sample}_R1.decontam.paired.fq.gz"
#    R2="$SRA_PATH/$sample/${sample}_R2.decontam.paired.fq.gz"
#    SAM="$ARG_DIR/mappings/${sample}_arg.sam"
#    BAM="$ARG_DIR/mappings/${sample}_arg.sorted.bam"
#    BAI="$BAM.bai"

    # Skip if both SAM and BAM index exist
#    if [ -f "$SAM" ] && [ -f "$BAI" ]; then
#      echo " ⚠ Skipping $sample: both SAM and BAM index already exist"
#      continue
#    fi

#    step_start=$(date +%s)
#    bowtie2 -x "$ARG_DIR/indices/arg_index" -1 "$R1" -2 "$R2" \
#      -S "$SAM" -p "$THREADS" --quiet
#    step_end=$(date +%s)

#    log_step "$sample" "Bowtie2_map" "$step_start" "$step_end" "$R1" "$SAM"
#    echo " ✔ Mapping complete for $sample"
#  done
#done


############### Convert SAM to BAM ####################################
conda activate samtools || exit 1
echo "[Step 4.5] Converting SAM to BAM"

find "$ARG_DIR/mappings" -type f -name "*.sam" | while read -r sam; do
  sample=$(basename "$sam" _arg.sam)
  bam="$ARG_DIR/mappings/${sample}_arg.sorted.bam"
  step_start=$(date +%s)

  samtools view -@ "$THREADS" -bS "$sam" | samtools sort -@ "$THREADS" -o "$bam"
  samtools index "$bam"

  step_end=$(date +%s)
  log_step "$sample" "SAM_to_BAM" "$step_start" "$step_end" "$sam" "$bam"
  echo "    ARG mapping and BAM complete for $sample"
done


############################################
# Step 5: Run CoverM  # Run this after you get all the sam to bam
############################################
#echo "[Step 5] Quantifying ARG abundance with CoverM"
#conda activate coverm || exit 1
#for bam in "$ARG_DIR/mappings"/*_arg.sorted.bam; do
#  sample=$(basename "$bam" _arg.sorted.bam)
#  step_start=$(date +%s)
#  coverm contig \
#    --bam-files "$bam" \
#    --methods rpkm tpm count \
#    --min-read-percent-identity 80 \
#    --min-read-aligned-percent 80 \
#    --threads "$THREADS" \
#    --output-file "$ARG_DIR/abundance/${sample}_arg_abundance.tsv"
#  step_end=$(date +%s)
#  log_step "$sample" "CoverM_ARG_abundance" "$step_start" "$step_end" "$bam" "$ARG_DIR/abundance/${sample}_arg_abundance.tsv"
#  echo "    Abundance profiling complete for $sample"
#done

#echo "✅ Pipeline complete: Global ARG abundance quantification ready."
