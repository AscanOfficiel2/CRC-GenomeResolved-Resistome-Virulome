#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=35:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/slurm-%j.out2
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/slurm-%j.err2

eval "$(conda shell.bash hook)"

project=$1
SRA_path="/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/DATA/$project"
collect="/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/DATA/collect_AMR_Env/$project"

fastq_screen_Before1="$collect/FQS_Before1"
fastq_screen_Before2="$collect/FQS_Before2"
quast_result="$collect/quast"
checkm2_result="$collect/checkm2"
Process_screen="$collect/Process"
Process_screen2="$collect/Process2"
NP="$collect/Nonpareil_$project"

mkdir -p "$fastq_screen_Before1" "$fastq_screen_Before2" "$quast_result" "$checkm2_result" "$Process_screen" "$Process_screen2" "$NP"

while read -r ids; do
    id_dir="$SRA_path/$ids"

    # FASTQ screen file collection
    r1_file="$id_dir/fastqscreen_${ids}/${ids}_R1.qc.nophix_screen.txt"
    r2_file="$id_dir/fastqscreen_${ids}/${ids}_R2.qc.nophix_screen.txt"
    [ -f "$r1_file" ] && cp "$r1_file" "$fastq_screen_Before1/"
    [ -f "$r2_file" ] && cp "$r2_file" "$fastq_screen_Before2/"

    # CheckM2 collection
    checkm_src="$id_dir/checkm2_${ids}/quality_report.tsv"
    checkm_dest="$checkm2_result/${ids}_quality_report.tsv"
    [ -f "$checkm_src" ] && cp "$checkm_src" "$checkm_dest"

    # QUAST collection
    quast_filtered="$id_dir/quast/quast_FILTERED/transposed_report.tsv"
    quast_simplified="$id_dir/quast/quast_SIMPLIFIED/transposed_report.tsv"
    [ -f "$quast_filtered" ] && cp "$quast_filtered" "$quast_result/${ids}_FILTERED_transposed_report.tsv"
    [ -f "$quast_simplified" ] && cp "$quast_simplified" "$quast_result/${ids}_SIMPLIFIED_transposed_report.tsv"

    # Nonpareil collection
    nonpareil_dir="$id_dir/nonpareil_${ids}"
    nonpareil_file="${nonpareil_dir}/${ids}_nonpareil.npo"
    nonpareil_dest="${NP}/${ids}_nonpareil.npo"
    if [ -f "$nonpareil_file" ]; then
        cp "$nonpareil_file" "$nonpareil_dest"
    fi

done < "$SRA_path/SRA_ids"

# Table processing for Fastq_screen R1
cd "$fastq_screen_Before1"
ls *_R1.qc.nophix_screen.txt > ids2.txt
for fname in $(cat ids2.txt); do
    SRR=$(basename "$fname" "_R1.qc.nophix_screen.txt")
    awk -F "\t" -v org="$SRR" '{print org"\t"$1"\t"$2"\t"$3}' "$fastq_screen_Before1/$SRR"_R1.qc.nophix_screen.txt > "$Process_screen/${SRR}_screen1.txt"
    sed -i '1d' "$Process_screen/${SRR}_screen1.txt"
    sed -i '1d' "$Process_screen/${SRR}_screen1.txt"
    sed -i 's/%//g' "$Process_screen/${SRR}_screen1.txt"
done

cd "$Process_screen"
for screen_file in *_screen1.txt; do
    SRR=$(basename "$screen_file" "_screen1.txt")
    while IFS= read -r ids; do
        printf "%s\t" "$ids" >> "$Process_screen/${SRR}_FQS_before1.txt"
    done < "$screen_file"
done
sed -i -e '$a\' $Process_screen/*_FQS_before1.txt
cat $Process_screen/*_FQS_before1.txt > "$fastq_screen_Before1/Fastqscreen_before1.csv"

# Table processing for Fastq_screen R2
cd "$fastq_screen_Before2"
ls *_R2.qc.nophix_screen.txt > ids2.txt
for fname in $(cat ids2.txt); do
    SRR=$(basename "$fname" "_R2.qc.nophix_screen.txt")
    awk -F "\t" -v org="$SRR" '{print org"\t"$1"\t"$2"\t"$3}' "$fastq_screen_Before2/$SRR"_R2.qc.nophix_screen.txt > "$Process_screen2/${SRR}_screen2.txt"
    sed -i '1d' "$Process_screen2/${SRR}_screen2.txt"
    sed -i '1d' "$Process_screen2/${SRR}_screen2.txt"
    sed -i 's/%//g' "$Process_screen2/${SRR}_screen2.txt"
done

cd "$Process_screen2"
for screen_file in *_screen2.txt; do
    SRR=$(basename "$screen_file" "_screen2.txt")
    while IFS= read -r ids; do
        printf "%s\t" "$ids" >> "$Process_screen2/${SRR}_FQS_before2.txt"
    done < "$screen_file"
done
sed -i -e '$a\' $Process_screen2/*_FQS_before2.txt
cat $Process_screen2/*_FQS_before2.txt > "$fastq_screen_Before2/Fastqscreen_before2.csv"

# MultiQC (optional)
multiqc "$fastq_screen_Before1" -o "$fastq_screen_Before1"
multiqc "$fastq_screen_Before2" -o "$fastq_screen_Before2"
