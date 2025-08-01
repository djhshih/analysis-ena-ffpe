#!/bin/bash

set -euo pipefail

root_outdir="../../../evaluations/depth/mpileup"
ref="../../../data/ref/Homo_sapiens_assembly38.fasta"

mkdir -p $root_outdir

for bam in ../../../data/bam/*.bam; do
    
    basename=$(basename $bam)
    filename=${basename%%.*}

    echo "Processing $basename..."

    # Include the flag below if the reference nucleotide needs to be labeled at that coordinate
    # This adds additional computational time and makes read base encoding more cryptic
    # Visit https://www.htslib.org/doc/samtools-mpileup.html for details
    # Flag: -f $ref

    samtools mpileup -d 0 $bam -o ${root_outdir}/${filename}_mpileup.tsv

    echo -e "\t pileup saved to ${root_outdir}/${filename}_mpileup.tsv \n"

done

echo "Done."