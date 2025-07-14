#!/bin/bash

njobs=72

indir=inputs
workflow=../../wdl/fastq_align_paired.wdl 
outdir=../../bam

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

if (( $npending > 0 )); then

	for input in $(head -n $njobs pending.vtr); do
		echo $input
		cromwell submit -i $indir/${input}.inputs $workflow

	done

fi
