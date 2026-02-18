#!/bin/bash
#SBATCH --job-name=assembly_binning
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=36:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/%x-%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/%x-%j.err

eval "$(conda shell.bash hook)" || true

# Inputs
project=${1:?Provide project accession, e.g. PRJEB6070}
SRA_path=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/DATA/$project
bbmap_path=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/bbmap

THREADS="${SLURM_CPUS_PER_TASK:-${SLURM_NTASKS:-56}}"
MEM_GB="${MEM_GB:-96}"

RUNTIME_LOG="$SRA_path/runtime_metrics_${project}_assembly_binning1.csv"
if [[ ! -f "$RUNTIME_LOG" ]]; then
  echo "Sample,Step,StartTime,EndTime,Duration_sec,Threads,Mem_GB,InputSize,OutputSize" > "$RUNTIME_LOG"
fi

cd "$SRA_path"
[[ -s "$SRA_path/SRA_ids1" ]] || { echo "Missing $SRA_path/SRA_ids"; exit 1; }

for ids in $(cat $SRA_path/SRA_ids1); do
  SAMPLE_DIR="$SRA_path/$ids"
  mkdir -p "$SAMPLE_DIR/assembly_qc" "$SAMPLE_DIR/quast" "$SAMPLE_DIR/depth" "$SAMPLE_DIR/bins"
  cd "$SAMPLE_DIR" || exit 1
  echo "== Processing $project / $ids =="

  log_step () {
    local sample=$1 step=$2 start=$3 end=$4 infile=$5 outfile=$6
    local dur=$((end-start))
    local insize=$( [ -f "$infile" ] && du -sh "$infile" | cut -f1 || echo "NA" )
    local outsize=$( [ -f "$outfile" ] && du -sh "$outfile" | cut -f1 || echo "NA" )
    echo "$sample,$step,$start,$end,$dur,$THREADS,$MEM_GB,$insize,$outsize" >> "$RUNTIME_LOG"
  }

  # Step 1: MEGAHIT
  R1="$SAMPLE_DIR/${ids}_R1.decontam.paired.fq.gz"
  R2="$SAMPLE_DIR/${ids}_R2.decontam.paired.fq.gz"
  SE="$SAMPLE_DIR/${ids}_singletons.decontam.fq.gz"
  ASSEMBLY_OUT="$SAMPLE_DIR/megahit"
  conda activate megahit
  step_start=$(date +%s)
  megahit --preset meta-large -1 "$R1" -2 "$R2" -r "$SE" -t "$THREADS" -o "$ASSEMBLY_OUT"
  step_end=$(date +%s)
  log_step "$ids" "MEGAHIT" "$step_start" "$step_end" "$R1" "$ASSEMBLY_OUT/final.contigs.fa"

  # Step 2: Simplify headers
  SIMPLIFIED="$SAMPLE_DIR/assembly_qc/${ids}_contigs_mod_headers.fasta"
  HEADER_MAP="$SAMPLE_DIR/assembly_qc/${ids}_contigs_header.map"
  CONTIGS_IN="$ASSEMBLY_OUT/final.contigs.fa"
  step_start=$(date +%s)
  perl /srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/used_scripts/MAG_script/simplifyFastaHeaders.pl \
    "$CONTIGS_IN" xx "$SIMPLIFIED" "$HEADER_MAP"
  step_end=$(date +%s)
  log_step "$ids" "SimplifyHeaders" "$step_start" "$step_end" "$CONTIGS_IN" "$SIMPLIFIED"

  # Step 3: Filter contigs >=1.5kb
  conda activate seqtk
  FILTERED="$SAMPLE_DIR/assembly_qc/${ids}_contigs_mod_headers.min1500.fasta"
  step_start=$(date +%s)
  seqtk seq -L 1500 "$SIMPLIFIED" > "$FILTERED"
  step_end=$(date +%s)
  log_step "$ids" "FilterContigs_>=1.5kb" "$step_start" "$step_end" "$SIMPLIFIED" "$FILTERED"

  # Step 4: QUAST
  conda activate quast
  QUAST_OUT1="$SAMPLE_DIR/quast/quast_SIMPLIFIED"
  QUAST_OUT2="$SAMPLE_DIR/quast/quast_FILTERED"
  step_start=$(date +%s)
  quast.py "$SIMPLIFIED" --threads "$THREADS" -o "$QUAST_OUT1"
  step_end=$(date +%s)
  log_step "$ids" "QUAST_SIMPLIFIED" "$step_start" "$step_end" "$SIMPLIFIED" "$QUAST_OUT1/report.txt"
  step_start=$(date +%s)
  quast.py "$FILTERED" --threads "$THREADS" -o "$QUAST_OUT2"
  step_end=$(date +%s)
  log_step "$ids" "QUAST_FILTERED" "$step_start" "$step_end" "$FILTERED" "$QUAST_OUT2/report.txt"

  # Step 5: Mapping & Binning
  ASSEMBLY="$FILTERED"
  step_start=$(date +%s)
  bash "$bbmap_path/bbmap.sh" ref="$ASSEMBLY" in1="$R1" in2="$R2" out="$SAMPLE_DIR/depth/${ids}_mapped.sam" threads=$THREADS
  conda activate samtools
  samtools view -bS "$SAMPLE_DIR/depth/${ids}_mapped.sam" | \
    samtools sort -@ "$THREADS" -o "$SAMPLE_DIR/depth/${ids}.sorted.bam"
  samtools index "$SAMPLE_DIR/depth/${ids}.sorted.bam"
  step_end=$(date +%s)
  log_step "$ids" "BBMap_Mapping" "$step_start" "$step_end" "$ASSEMBLY" "$SAMPLE_DIR/depth/${ids}.sorted.bam"

  conda activate metabat2
  DEPTH_FILE="$SAMPLE_DIR/depth/${ids}.depth.txt"
  step_start=$(date +%s)
  jgi_summarize_bam_contig_depths --outputDepth "$DEPTH_FILE" "$SAMPLE_DIR/depth/${ids}.sorted.bam"
  step_end=$(date +%s)
  log_step "$ids" "Compute_Depth" "$step_start" "$step_end" "$SAMPLE_DIR/depth/${ids}.sorted.bam" "$DEPTH_FILE"

  BIN_DIR="$SAMPLE_DIR/bins"
  step_start=$(date +%s)
  metabat2 -i "$ASSEMBLY" -a "$DEPTH_FILE" -o "$BIN_DIR/${ids}_bin" -t "$THREADS" --minContig 1500 --maxP 90 --minS 80 --maxEdges 200 --minCV 1
  step_end=$(date +%s)
  log_step "$ids" "MetaBAT2" "$step_start" "$step_end" "$ASSEMBLY" "$BIN_DIR"

  echo "[$ids] Assembly + Binning complete"


  ###############################
## CheckM2 on bins
###############################
conda activate checkm2
CHECKM_OUT="$SAMPLE_DIR/checkm2_$ids"
mkdir -p "$CHECKM_OUT"
export CHECKM2DB="/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Shotgun-metagenomics/checkm2_databases/CheckM2_database/"
step_start=$(date +%s)
checkm2 predict -x fa -i "$BIN_DIR" -o "$CHECKM_OUT" --threads "$THREADS" --database_path "$CHECKM2DB/uniref100.KO.1.dmnd"
step_end=$(date +%s)
log_step "$ids" "CheckM2" "$step_start" "$step_end" "$BIN_DIR" "$CHECKM_OUT"


done

