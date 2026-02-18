#!/usr/bin/env bash
## Last update 23-08-2025 by Suleiman & AbdulAziz
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=35:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/slurm-%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/slurm-%j.err



eval "$(conda shell.bash hook)" || true

# ------------------ inputs & paths ------------------
project=${1:?Provide project accession, e.g. SRP099122}

sh_path=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/used_scripts/MAG_script
SRA_path=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/DATA/$project
bbmap_path=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/bbmap
adapters=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/bbmap/resources
export FASTQ_SCREEN_CONF="/home/suleiman.aminu/anaconda3/envs/fastq_screen/share/fastq-screen-0.15.3-0/fastq_screen.conf"

[[ -f "$FASTQ_SCREEN_CONF" ]] || { echo "FASTQ_SCREEN_CONF not found: $FASTQ_SCREEN_CONF" >&2; exit 1; }

THREADS="${SLURM_CPUS_PER_TASK:-${SLURM_NTASKS:-56}}"
MEM_GB="${MEM_GB:-96}"
MINLEN="${MINLEN:-50}"

# CSV runtime log file
RUNTIME_LOG="$SRA_path/runtime_metrics_${project}.csv"
echo "Sample,Step,StartTime,EndTime,Duration_sec,Threads,Mem_GB,InputSize,OutputSize" > "$RUNTIME_LOG"

# ------------------ start ------------------
cd "$SRA_path"
[[ -s "$SRA_path/SRA_ids" ]] || { echo "Missing $SRA_path/SRA_ids"; exit 1; }

# Add SRA toolkit to PATH
export PATH=$PATH:/home/suleiman.aminu/sratoolkit.3.1.0-ubuntu64/bin

# ------------------ Loop over all samples ------------------
for ids in $(cat $SRA_path/SRA_ids); do
  SAMPLE_DIR="$SRA_path/$ids"
  mkdir -p "$SAMPLE_DIR"
  cd "$SAMPLE_DIR" || exit 1
  echo "== Processing $project / $ids =="

  ######################################
  # Function to log step runtimes
  ######################################
  log_step () {
    local sample=$1
    local step=$2
    local start=$3
    local end=$4
    local infile=$5
    local outfile=$6
    local dur=$((end-start))
    local insize=$( [ -f "$infile" ] && du -sh "$infile" | cut -f1 || echo "NA" )
    local outsize=$( [ -f "$outfile" ] && du -sh "$outfile" | cut -f1 || echo "NA" )
    echo "$sample,$step,$start,$end,$dur,$THREADS,$MEM_GB,$insize,$outsize" >> "$RUNTIME_LOG"
  }

  # -------------------- 1. Fetch FASTQ --------------------
  step_start=$(date +%s)
  if [[ ! -f "${ids}_1.fastq.gz" || ! -f "${ids}_2.fastq.gz" ]]; then
    prefetch "$ids" -O .
    fasterq-dump "$ids" --split-files --threads "$THREADS" -O .
    pigz -p "$THREADS" "${ids}_1.fastq" "${ids}_2.fastq"
  fi
  step_end=$(date +%s)
  log_step "$ids" "Fetch_FASTQ" $step_start $step_end "${ids}.sra" "${ids}_1.fastq.gz"

  # -------------------- 2. Adapter trimming --------------------
  R1_QC="${ids}_R1.qc.fq.gz"
  R2_QC="${ids}_R2.qc.fq.gz"
  step_start=$(date +%s)
  if [[ ! -f "$R1_QC" || ! -f "$R2_QC" ]]; then
    bash "$bbmap_path/bbduk.sh" in1="${ids}_1.fastq.gz" in2="${ids}_2.fastq.gz" \
      out1="$R1_QC" out2="$R2_QC" ref="${adapters}/adapters.fa" \
      ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=20 minlen="$MINLEN" threads="$THREADS"
  fi
  step_end=$(date +%s)
  log_step "$ids" "AdapterTrim" $step_start $step_end "${ids}_1.fastq.gz" "$R1_QC"

  # -------------------- 3. PhiX removal --------------------
  R1_NOPHIX="${ids}_R1.qc.nophix.fq.gz"
  R2_NOPHIX="${ids}_R2.qc.nophix.fq.gz"
  step_start=$(date +%s)
  if [[ ! -f "$R1_NOPHIX" || ! -f "$R2_NOPHIX" ]]; then
    bash "$bbmap_path/bbduk.sh" in1="$R1_QC" in2="$R2_QC" \
      out1="$R1_NOPHIX" out2="$R2_NOPHIX" ref="${adapters}/phix174_ill.ref.fa.gz" \
      k=31 hdist=1 stats="bbduk_phix_${ids}.txt" threads="$THREADS"
  fi
  step_end=$(date +%s)
  log_step "$ids" "PhiXRemoval" $step_start $step_end "$R1_QC" "$R1_NOPHIX"

  # -------------------- 4. FastQ Screen (host check) --------------------
  conda activate fastq_screen
  outd="$SAMPLE_DIR/fastqscreen_${ids}"
  mkdir -p "$outd"
  step_start=$(date +%s)
  fastq_screen --aligner bowtie --threads "$THREADS" --subset 0 \
    --conf "$FASTQ_SCREEN_CONF" --outdir "$outd" --tag "$R1_NOPHIX"
  fastq_screen --aligner bowtie --threads "$THREADS" --subset 0 \
    --conf "$FASTQ_SCREEN_CONF" --outdir "$outd" --tag "$R2_NOPHIX"
  step_end=$(date +%s)
  log_step "$ids" "FastQScreen" $step_start $step_end "$R1_NOPHIX" "$outd/${ids}_screen.txt"

  # -------------------- 5. Repair reads --------------------

  # Strip .fq.gz → add .tagged.fastq.gz
  R1_TAG="$outd/$(basename "${R1_NOPHIX%.fq.gz}").tagged.fastq.gz"
  R2_TAG="$outd/$(basename "${R2_NOPHIX%.fq.gz}").tagged.fastq.gz"

  # Keep only reads with no hits
  fastq_screen --nohits "$R1_TAG"
  fastq_screen --nohits "$R2_TAG"

  # Locate the filtered FastQ files in fastq_screen output dir
  R1_FILT="$SAMPLE_DIR/$(basename "${R1_NOPHIX%.fq.gz}").tagged.tagged_filter.fastq.gz"
  R2_FILT="$SAMPLE_DIR/$(basename "${R2_NOPHIX%.fq.gz}").tagged.tagged_filter.fastq.gz"

  # Move them into SAMPLE_DIR for consistency
  cp "$R1_FILT" "$SAMPLE_DIR/"
  cp "$R2_FILT" "$SAMPLE_DIR/"


  # If gzipped, fix names
  #[[ -f "${R1_FILT}.gz" ]] && R1_FILT="${R1_FILT}.gz"
  #[[ -f "${R2_FILT}.gz" ]] && R2_FILT="${R2_FILT}.gz"

  # Repair pairs and gzip
 step_start=$(date +%s)
  bash "$bbmap_path/repair.sh" in1="$R1_FILT" in2="$R2_FILT" \
    out1="$SAMPLE_DIR/${ids}_R1.decontam.paired.fq.gz" \
    out2="$SAMPLE_DIR/${ids}_R2.decontam.paired.fq.gz" \
    outs="$SAMPLE_DIR/${ids}_singletons.decontam.fq.gz" \
    overwrite=true
  step_end=$(date +%s)
  log_step "$ids" "Repair" $step_start $step_end "$R1_FILT" "$SAMPLE_DIR/${ids}_R1.decontam.paired.fq.gz"

  # -------------------- 6. FastQC --------------------
  conda activate fastqc
  PRE_DIR="$SAMPLE_DIR/fastqc_prehost"
  POST_DIR="$SAMPLE_DIR/fastqc_posthost"
  mkdir -p "$PRE_DIR" "$POST_DIR"
  step_start=$(date +%s)
  fastqc -t "$THREADS" -o "$PRE_DIR" "$R1_NOPHIX" "$R2_NOPHIX"
  fastqc -t "$THREADS" -o "$POST_DIR" \
    "$SAMPLE_DIR/${ids}_R1.decontam.paired.fq.gz" \
    "$SAMPLE_DIR/${ids}_R2.decontam.paired.fq.gz"
  step_end=$(date +%s)
  log_step "$ids" "FastQC" $step_start $step_end "$R1_NOPHIX" "$SAMPLE_DIR/${ids}_R1.decontam.paired.fq.gz"

  # -------------------- 7. Nonpareil coverage --------------------
  conda activate nonpareil
  NONP_OUT="$SAMPLE_DIR/nonpareil_${ids}"
  mkdir -p "$NONP_OUT"
  step_start=$(date +%s)
  nonpareil -s "$SAMPLE_DIR/${ids}_R1.decontam.paired.fq.gz" \
            -T kmer -f fastq \
            -b "$NONP_OUT/${ids}_nonpareil" \
            -t "$THREADS"
  step_end=$(date +%s)
  log_step "$ids" "Nonpareil" $step_start $step_end "$SAMPLE_DIR/${ids}_R1.decontam.paired.fq.gz" "$NONP_OUT/${ids}_nonpareil.nps"

done

