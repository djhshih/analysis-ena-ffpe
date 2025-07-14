#!/bin/bash

#PBS -V
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=4,mem=10g,walltime=08:00:00
#PBS -q cgsd
#PBS -o /home/wtjtse/projects/analysis-ccoc-ts/submit/fastq_align_paired/test/test.stdout
#PBS -e /home/wtjtse/projects/analysis-ccoc-ts/submit/fastq_align_paired/test/test.stderr
#PBS -N test


njobs=100

indir=$(realpath inputs)
workflow=$(realpath ../../wdl/fastq_align_paired.wdl)
outdir=$(realpath ../../bam)

mkdir -p log
mkdir -p job

logdir=$(realpath log)

function fstem {
	sed -E 's|.*/||;s|[.].*||'
}

# update hard links to output files
(cd $outdir && ./get.sh)

# find pending samples that still need to be processed
comm -2 -3 <(ls -1 $indir/*.inputs | fstem | sort) <(ls -1 $outdir/*.bam | fstem | sort) \
	> pending.vtr

# number of pending samples
npending=$(wc -l pending.vtr | cut -d ' ' -f 1)
echo "number of samples pending: $npending" >&2

cromwell_path=/home/wtjtse/local/opt/cromwell/bin/cromwell

if (( $npending > 0 )); then

	for input in $(head -n $njobs pending.vtr); do
		echo $input
		echo "
#PBS -V 
#PBS -S /bin/bash 
#PBS -l nodes=1:ppn=8,mem=10g,walltime=24:00:00
#PBS -q cgsd
#PBS -o ${logdir}/${input}.stdout
#PBS -e ${logdir}/${input}.stderr
#PBS -N ${input}
module load bwa 
module load samtools
module load java
module load gcc
module load sambamba" > job/${input}.sh
	
		echo ${cromwell_path} run -i $indir/${input}.inputs $workflow >> job/${input}.sh
		 
	done

fi
