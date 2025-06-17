######Make index. files for reference genome

mkdir -p bowtie2_index
time apptainer exec /work2/03302/lconcia/sif_files/bowtie2_2.4.4.sif \
bowtie2-build --threads 48 Zm-B73-REFERENCE-NAM-5.0_without_scaffold.fa bowtie2_index/Zm-B73-REFERENCE-NAM-5.0
