######Make index. files for reference genome**
#!/bin/bash
#SBATCH -J bowtie2-build_chromHMM.job # job name
#SBATCH -o bowtie2-build_chromHMM.%j.out # output file name (%j expands to jobID)
#SBATCH -e bowtie2-build_chromHMM.%j.err # error file name (%j expands to jobID)
#SBATCH -N 1 # number of nodes requested ## up to 96 for SKX
#SBATCH -c 48 # --cpus-per-task ## hicpro config (fw / rv)
#SBATCH --ntasks-per-node 1 # --ntasks-per-node
#SBATCH --mem=0 # used all the memory available on the node
#SBATCH -p spr # queue (partition)
#SBATCH -t 48:00:00 # run time (hh:mm:ss)
#SBATCH -A MCB24009 # Allocation name to charge job against
module load tacc-apptainer/1.2.5
module unload xalt
./bowtie2-build_chromHMM.sh
#### bowtie2-build_chromHMM.sh
mkdir -p bowtie2_index
time apptainer exec /work2/03302/lconcia/sif_files/bowtie2_2.4.4.sif \
bowtie2-build --threads 48 Zm-B73-REFERENCE-NAM-5.0_without_scaffold.fa
bowtie2_index/Zm-B73-REFERENCE-NAM-5.0
#####CHIP-seq read alignment
#!/bin/bash
#SBATCH -J CHIP-seq_alignment.job # job name
#SBATCH -o CHIP-seq_alignment.%j.out # output file name (%j expands to jobID)
#SBATCH -e CHIP-seq_alignment.%j.err # error file name (%j expands to jobID)
#SBATCH -N 1 # number of nodes requested ## up to 96 for SKX
#SBATCH -c 48 # --cpus-per-task ## hicpro config (fw / rv)
#SBATCH --ntasks-per-node 1 # --ntasks-per-node
#SBATCH --mem=0 # used all the memory available on the node
#SBATCH -p spr # queue (partition)
#SBATCH -t 48:00:00 # run time (hh:mm:ss)
#SBATCH -A MCB24009 # Allocation name to charge job against
module load tacc-apptainer/1.2.5
module unload xalt
# Create a clean output directory
mkdir -p aligned_bams
# Loop through all _1.fastq files only
for r1 in *_1.fastq; do
r2=${r1/_1.fastq/_2.fastq}
sample=$(basename $r1 _1.fastq)
echo "Processing $sample ..."
# Run alignment and generate BAM

time apptainer exec /work2/03302/lconcia/sif_files/bowtie2_2.4.4.sif bowtie2 -p 48 --very-
sensitive \

-x bowtie2_index/Zm-B73-REFERENCE-NAM-5.0 \
-1 $r1 -2 $r2 | \

apptainer exec /work2/03302/lconcia/sif_files/samtools_v1.9-4-deb_cv1.sif samtools view -@
48 -bS - > aligned_bams/${sample}_prefiltered.bam
done
#### remove duplicates
#!/bin/bash
#SBATCH -J rmduplicates_chipseq.job # job name
#SBATCH -o rmduplicates_chipseq.%j.out # output file name (%j expands to jobID)
#SBATCH -e rmduplicates_chipseq.%j.err # error file name (%j expands to jobID)
#SBATCH -N 1 # number of nodes requested ## up to 96 for SKX
#SBATCH -c 48 # --cpus-per-task ## hicpro config (fw / rv)
#SBATCH --ntasks-per-node 1 # --ntasks-per-node
#SBATCH --mem=0 # used all the memory available on the node
#SBATCH -p spr # queue (partition)
#SBATCH -t 48:00:00 # run time (hh:mm:ss)
#SBATCH -A MCB24009 # Allocation name to charge job against
module load tacc-apptainer/1.2.5
module unload xalt
mkdir -p dedup_bams
for bam in aligned_bams/*_prefiltered.bam; do
sample=$(basename $bam _prefiltered.bam)
echo "Processing $sample ..."
# Step 1: Sort by read name (needed for fixmate)
apptainer exec /work2/03302/lconcia/sif_files/samtools_v1.9-4-deb_cv1.sif \
samtools sort -n -@ 48 -o dedup_bams/${sample}_namesorted.bam $bam
# Step 2: Run fixmate to add MC tags
apptainer exec /work2/03302/lconcia/sif_files/samtools_v1.9-4-deb_cv1.sif \
samtools fixmate -m dedup_bams/${sample}_namesorted.bam
dedup_bams/${sample}_fixmate.bam
# Step 3: Sort by coordinate
apptainer exec /work2/03302/lconcia/sif_files/samtools_v1.9-4-deb_cv1.sif \
samtools sort -@ 48 -o dedup_bams/${sample}_sorted.bam
dedup_bams/${sample}_fixmate.bam
# Step 4: Mark and remove duplicates
apptainer exec /work2/03302/lconcia/sif_files/samtools_v1.9-4-deb_cv1.sif \
samtools markdup -r -@ 48 dedup_bams/${sample}_sorted.bam
dedup_bams/${sample}_dedup.bam
# Step 5: Index final deduplicated BAM
apptainer exec /work2/03302/lconcia/sif_files/samtools_v1.9-4-deb_cv1.sif \
samtools index dedup_bams/${sample}_dedup.bam
done
#### bamtobed_chromHMM
#!/bin/bash
#SBATCH -J bamtobed_chromHMM.job # job name
