#!/bin/bash
#SBATCH --job-name=wdl-variant-call-mutect2       # 1. Job name
#SBATCH --mail-type=BEGIN,END,FAIL    # 2. Send email upon events (Options: NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=moyukh@hku.hk     #    Email address to receive notification
#SBATCH --partition=amd               # 3. Request a partition
#SBATCH --qos=normal                  # 4. Request a QoS
#SBATCH --ntasks=80                   # 5. Request total number of tasks (MPI workers)
#SBATCH --nodes=1                     #    Request number of node(s)
#SBATCH --mem=256G                     # 6. Request total amount of RAM
#SBATCH --time=2-00:00:00            # 7. Job execution duration limit day-hour:min:sec
#SBATCH --output=sbatch_log/%x_%j.out            # 8. Standard output log as $job_name_$job_id.out
#SBATCH --error=sbatch_log/%x_%j.err             #    Standard error log as $job_name_$job_id.err

# print the start time
date

echo -e "\nSLURM_NTASKS: $SLURM_NTASKS\n"

module load gnuparallel

# switch to proper environment
source ~/.bashrc
conda activate ena-ffpe

ls ../../data/vcf_filtermutectcalls_obf/*/*.vcf > dlazy_samples.txt
djobs dlazy_samples.txt bash filter_alignment_artifacts.sh

pdlazy job -j $(($SLURM_NTASKS/2))

# print the end time
date
