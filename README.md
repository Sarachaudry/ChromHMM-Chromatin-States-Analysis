######Make index. files for reference genome

#!/bin/bash
#SBATCH -J  bowtie2-build_chromHMM.job          # job name
#SBATCH -o  bowtie2-build_chromHMM.%j.out  # output file name (%j expands to jobID)
#SBATCH -e  bowtie2-build_chromHMM.%j.err  # error  file name (%j expands to jobID)
#SBATCH -N 1                                                # number of nodes requested ## up to 96 for SKX
#SBATCH -c 48                                               # --cpus-per-task   ## hicpro config (fw / rv)
#SBATCH --ntasks-per-node 1                                 # --ntasks-per-node
#SBATCH --mem=0                                             # used all the memory available on the node
#SBATCH -p spr                                              # queue (partition)
#SBATCH -t 48:00:00                                         # run time (hh:mm:ss)
#SBATCH -A MCB24009                                        # Allocation name to charge job against
module load tacc-apptainer/1.2.5
module unload xalt
./bowtie2-build_chromHMM.sh
