######Make index. files for reference genome

mkdir -p bowtie2_index
time apptainer exec /work2/03302/lconcia/sif_files/bowtie2_2.4.4.sif \
bowtie2-build --threads 48 Zm-B73-REFERENCE-NAM-5.0_without_scaffold.fa bowtie2_index/Zm-B73-REFERENCE-NAM-5.0

#####CHIP-seq read alignment

# Create a clean output directory
mkdir -p aligned_bams

# Loop through all _1.fastq files only
for r1 in *_1.fastq; do
    r2=${r1/_1.fastq/_2.fastq}
    sample=$(basename $r1 _1.fastq)

    echo "Processing $sample ..."

    # Run alignment and generate BAM
    time apptainer exec /work2/03302/lconcia/sif_files/bowtie2_2.4.4.sif bowtie2 -p 48 --very-sensitive \
        -x bowtie2_index/Zm-B73-REFERENCE-NAM-5.0 \
        -1 $r1 -2 $r2 | \
    apptainer exec /work2/03302/lconcia/sif_files/samtools_v1.9-4-deb_cv1.sif samtools view -@ 48 -bS - > aligned_bams/${sample}_prefiltered.bam
done

#### remove duplicates

# Create a clean output directory
mkdir -p dedup_bams

for bam in aligned_bams/*_prefiltered.bam; do
    sample=$(basename $bam _prefiltered.bam)

    echo "Processing $sample ..."

    # Step 1: Sort by read name (needed for fixmate)
    apptainer exec /work2/03302/lconcia/sif_files/samtools_v1.9-4-deb_cv1.sif \
    samtools sort -n -@ 48 -o dedup_bams/${sample}_namesorted.bam $bam

    # Step 2: Run fixmate to add MC tags
    apptainer exec /work2/03302/lconcia/sif_files/samtools_v1.9-4-deb_cv1.sif \
    samtools fixmate -m dedup_bams/${sample}_namesorted.bam dedup_bams/${sample}_fixmate.bam

    # Step 3: Sort by coordinate
    apptainer exec /work2/03302/lconcia/sif_files/samtools_v1.9-4-deb_cv1.sif \
    samtools sort -@ 48 -o dedup_bams/${sample}_sorted.bam dedup_bams/${sample}_fixmate.bam

    # Step 4: Mark and remove duplicates
    apptainer exec /work2/03302/lconcia/sif_files/samtools_v1.9-4-deb_cv1.sif \
    samtools markdup -r -@ 48 dedup_bams/${sample}_sorted.bam dedup_bams/${sample}_dedup.bam

    # Step 5: Index final deduplicated BAM
    apptainer exec /work2/03302/lconcia/sif_files/samtools_v1.9-4-deb_cv1.sif \
    samtools index dedup_bams/${sample}_dedup.bam

done

#### bamtobed_chromHMM

for bam in dedup_bams/*_dedup.bam; do
    sample=$(basename $bam _dedup.bam)

    echo "Filtering and converting $sample to BED..."

    apptainer exec /work2/03302/lconcia/sif_files/samtools_v1.9-4-deb_cv1.sif samtools view -@ 48 -b -F 4 -q 30 $bam | \
    apptainer exec /work2/03302/lconcia/sif_files/bedtools_2.31.0.sif bedtools bamtobed -i - > bedfiles_for_chromHMM/${sample}_filtered.bed
done

#### chromHMM job

# Run ChromHMM inside container using correct internal path
apptainer exec /work2/03302/lconcia/sif_files/chromhmm_1.26--hdfd78af_0.sif \
java -mx4000M -jar /usr/local/share/chromhmm-1.26-0/ChromHMM.jar BinarizeBed \
  -b 200 \
  /work2/03302/lconcia/references/maize/assemblies/B73/Zm-B73-REFERENCE-NAM-5.0.sizes.genome \
  /scratch/99999/ha20be/Akram_Basslab_PhD_HiC_analyses_revised/chromHMM_revised_manuscript/sra_download/bedfiles_for_chromHMM \
  cellmarkfiletable.txt \
  Binarized_RootTip

#### chromHMM learn model job

apptainer exec /work2/03302/lconcia/sif_files/chromhmm_1.26--hdfd78af_0.sif \
java -mx8000M -jar /usr/local/share/chromhmm-1.26-0/ChromHMM.jar LearnModel \
  Binarized_RootTip \
  ChromHMM_RootTip_10state \
  10 \
  B73v5

**####change marks order on x-axis for emission plots**

apptainer exec /work2/03302/lconcia/sif_files/chromhmm_1.26--hdfd78af_0.sif \
java -mx4000M -jar /usr/local/share/chromhmm-1.26-0/ChromHMM.jar \
Reorder -f emissionPlot_order -stateordering emission \
ChromHMM_RootTip_9state/model_9.txt \
ChromHMM_RootTip_9state/Reordered_9state_both

####loop

#!/bin/bash

SIF_PATH="/work2/03302/lconcia/sif_files/chromhmm_1.26--hdfd78af_0.sif"
ORDER_FILE="emissionPlot_order"

# Loop over state directories
for NSTATE in 3 4 5 
do
    MODEL_DIR="ChromHMM_RootTip_${NSTATE}state"
    MODEL_FILE="${MODEL_DIR}/model_${NSTATE}.txt"
    OUT_DIR="${MODEL_DIR}/Reordered_${NSTATE}state_both"

    # Create output directory if not exists
    mkdir -p "$OUT_DIR"

    # Run ChromHMM Reorder
    apptainer exec "$SIF_PATH" \
    java -mx4000M -jar /usr/local/share/chromhmm-1.26-0/ChromHMM.jar \
    Reorder \
    -f "$ORDER_FILE" \
    -stateordering emission \
    "$MODEL_FILE" \
    "$OUT_DIR"
done
