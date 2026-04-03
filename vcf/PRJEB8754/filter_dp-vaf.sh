#!/bin/bash

set -euo pipefail

root_indir="dup-unmarked_filtered_pass-orient-pos-sb"
root_outdir="dup-unmarked_filtered_pass-orient-pos-sb-vaf-dp"
mkdir -p $root_outdir

filter_expression='(FMT/AD[0:0] + FMT/AD[0:1]) >= 100 && (FMT/AD[0:1] / (FMT/AD[0:0] + FMT/AD[0:1])) >= 0.02'
echo -e "Filtering Expression: $filter_expression"

for vcf in $root_indir/*/*.vcf.gz; do
    
    filename=$(basename $vcf)
    sample_name=${filename%%.*}

    echo $sample_name

    outdir="${root_outdir}/${sample_name}"
    mkdir -p $outdir

    echo -e "\nFiltering $filename"
    bcftools view -i "$filter_expression" $vcf -Oz -o "${outdir}/${sample_name}.vcf.gz"
    echo "Filtered VCF saved to ${outdir}/${sample_name}.vcf.gz"

    echo "Indexing ${outdir}/${sample_name}.vcf.gz"
    bcftools index -t "${outdir}/${sample_name}.vcf.gz"

done

echo -e "\nFinished."
