### 1) Make index. files for reference genome
```bash
mkdir -p bowtie2_index
bowtie2-build --threads 48 Zm-B73-REFERENCE-NAM-5.0_without_scaffold.fa \
bowtie2_index/Zm-B73-REFERENCE-NAM-5.0
```
### 2) ChIP-seq read alignment

##### Create the output directory for aligned bam files
```bash
mkdir -p aligned_bams
```
##### Loop through all \*_1.fastq files files

```bash
for r1 in *_1.fastq;
do
    r2=${r1/_1.fastq/_2.fastq}
    sample=$(basename $r1 _1.fastq)
    
    echo "Processing $sample ..."
    
    # Run alignment and generate BAM
    bowtie2 -p 48 --very-sensitive \
    -x bowtie2_index/Zm-B73-REFERENCE-NAM-5.0 \
    -1 $r1 -2 $r2 | \
    samtools view -@ 48 -bS - > aligned_bams/${sample}_aligned.bam
done
```

### 3) Remove duplicated reads

##### Create the output directory for de-duplicated bam files
```bash
mkdir -p dedup_bams
```
##### Loop through all \*_prefiltered.bam files

```bash
for bam in aligned_bams/*_aligned.bam;
do
    sample=$(basename $bam _aligned.bam)

    echo "Processing $sample ..."

    # Step 1: Sort by read name (needed for fixmate)
    samtools sort -n -@ 48 -o dedup_bams/${sample}_namesorted.bam $bam

    # Step 2: Run fixmate to add MC tags
    samtools fixmate -m dedup_bams/${sample}_namesorted.bam dedup_bams/${sample}_fixmate.bam

    # Step 3: Sort by coordinate
    samtools sort -@ 48 -o dedup_bams/${sample}_fixmate.bam dedup_bams/${sample}_sorted.bam

    # Step 4: Mark and remove duplicates
    samtools markdup -r -@ 48 dedup_bams/${sample}_sorted.bam dedup_bams/${sample}_dedup.bam

    # Step 5: Index final deduplicated BAM
    samtools index dedup_bams/${sample}_dedup.bam

done
```

### 4) Convert to bam files to bed for chromHMM
##### Create the output directory for bed files
```bash
mkdir -p bedfiles_for_chromHMM
```
##### Loop through all \*_dedup.bam files 
```bash
for bam in dedup_bams/*_dedup.bam;
do
    sample=$(basename $bam _dedup.bam)

    echo "Filtering and converting $sample to BED..."

    samtools view -@ 48 -b -F 4 -q 30 $bam | \
    bedtools bamtobed -i - > bedfiles_for_chromHMM/${sample}_filtered.bed
done
```
### 5) Run chromHMM
##### Generate binarized data files from bed files
```bash
java -mx4000M -jar /usr/local/share/chromhmm-1.26-0/ChromHMM.jar BinarizeBed \
    -b 200 \
    /references/maize/assemblies/B73/Zm-B73-REFERENCE-NAM-5.0.sizes.genome \
    bedfiles_for_chromHMM \
   cellmarkfiletable.txt \
    Binarized_RootTip
```
##### Learn chromatin state models
```bash
java -mx8000M -jar /usr/local/share/chromhmm-1.26-0/ChromHMM.jar LearnModel \
    Binarized_RootTip \
    ChromHMM_RootTip_10state \
    10 \
    B73v5
```
##### Sort the marks on the x-axis in the emission plots on the results of one run
```bash
java -mx4000M -jar /usr/local/share/chromhmm-1.26-0/ChromHMM.jar \
Reorder -f emissionPlot_order -stateordering emission \
ChromHMM_RootTip_9state/model_9.txt \
ChromHMM_RootTip_9state/Reordered_9state_both
```
##### Repeat the sorting of marks on all the results
```bash
#!/bin/bash

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
    java -mx4000M -jar /usr/local/share/chromhmm-1.26-0/ChromHMM.jar \
    Reorder \
    -f "$ORDER_FILE" \
    -stateordering emission \
    "$MODEL_FILE" \
    "$OUT_DIR"
done
```
