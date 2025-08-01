#!/bin/bash

set -euo pipefail

root_outdir="../../data/vcf_pass-only"

mkdir -p $root_outdir

for vcf in ../../data/vcf_filtermutectcalls_obf/*/*.vcf; do
    
    filename=$(basename $vcf)
    sample_name=${filename%%.*}

    echo $sample_name

    outdir="${root_outdir}/${sample_name}"
    mkdir -p $outdir

    bcftools view -i 'FILTER="PASS"' -o "${outdir}/${sample_name}.vcf" $vcf

    # bcftools view -i 'FMT/AD[0:0] + FMT/AD[0:1] >= 10 && FMT/AD[0:1]'

done