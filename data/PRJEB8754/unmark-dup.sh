#!/usr/bin/env bash

set -euo pipefail

bam_indir="bam"
outdir_root="bam_dup-unmarked"

for bam in $bam_indir/*/*.bam; do

    echo Unmarking duplicates from: $bam

    sample=$(basename $(dirname $bam))
    outdir=$outdir_root/$sample
    mkdir -p $outdir

    gatk UnmarkDuplicates \
        -I $bam \
        -O $outdir/$sample.bam

    samtools index -@4 -b $outdir/$sample.bam -o $outdir/$sample.bai

done
