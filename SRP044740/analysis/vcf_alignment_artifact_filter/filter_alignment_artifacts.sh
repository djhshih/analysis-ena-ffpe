#!/bin/bash

set -euo pipefail

vcf=$1
if [[ -z $1 ]]; then
    echo "Usage: $0 <vcf_file>"
    exit 1
fi

# vcf_paths=../../data/vcf_filtermutectcalls_obf/*/*.vcf
referece_genome="../../data/ref/Homo_sapiens_assembly38.fasta"
bwa_mem_index_image="../../data/gatk-test-data/mutect2/Homo_sapiens_assembly38.index_bundle"

outdir_root="../../data/vcf_aaf_obf"

basename=$(basename $vcf)
sample_name=${basename%%.*}

echo "Processing $basename with GATK Alignment Artifact Filter"

bam="../../data/bam/${sample_name}.bam"
if [[ ! -f $bam ]]; then
    echo "BAM file does not exist at: $bam"
    exit 2
fi

outdir="${outdir_root}/${sample_name}"
mkdir -p $outdir

out_vcf="${outdir}/${sample_name}.vcf"

gatk FilterAlignmentArtifacts \
    -R $referece_genome \
    -V $vcf \
    -I $bam \
    --bwa-mem-index-image $bwa_mem_index_image \
    --dont-skip-filtered-variants \
    -O $out_vcf

echo -e "\n\tComplete."
