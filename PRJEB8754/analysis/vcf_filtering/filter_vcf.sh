#!/bin/env/bash

set -euo pipefail

vcf_dir="../../data/vcf_filtermutectcalls"
root_outdir="../../data/vcf_ad_filtered"

for vcf in ${vcf_dir}/*/*.vcf; do

    filename=$(basename $vcf)
    sample_id=${filename%%.*}
    
    echo $sample_id
    outdir=${root_outdir}/${sample_id}
    mkdir -p $outdir

    bcftools view -o ${outdir}/${sample_id}.vcf --include 'FMT/AD[0:0] + FMT/AD[0:1] >= 10 && FMT/AD[0:1] >= 3' $vcf

done