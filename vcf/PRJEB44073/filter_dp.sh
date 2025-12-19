#!/bin/bash

set -euo pipefail

root_outdir="vcf_filtered_pass-orientation-dp10"
mkdir -p $root_outdir

filter_expression='FMT/DP>=10'
echo -e "Filtering Expression: $filter_expression"

for vcf in vcf_filtered_pass_orientation/*/*.vcf; do
    
    filename=$(basename $vcf)
    sample_name=${filename%%.*}

    echo $sample_name
    
    outdir="${root_outdir}/${sample_name}"
    mkdir -p $outdir

    echo -e "\nFiltering $filename"
    bcftools view -i "$filter_expression" $vcf -o "${outdir}/${sample_name}.vcf"
    echo "Filtered VCF saved to ${outdir}/${sample_name}.vcf"

    echo "Indexing ${outdir}/${sample_name}.vcf"
    gatk IndexFeatureFile -I "${outdir}/${sample_name}.vcf"

done

echo -e "\nFinished."