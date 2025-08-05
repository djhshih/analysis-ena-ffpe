#!/bin/bash

set -euo pipefail

root_outdir="../../data/vcf_pass-orient-pos-sb_ad_filtered"
mkdir -p $root_outdir

filter_expression='(FILTER="PASS" | FILTER="strand_bias" | FILTER="orientation" | FILTER="position" | FILTER="strand_bias;orientation" | FILTER="strand_bias;position" | FILTER="orientation;position" | FILTER="strand_bias;orientation;position") & (FMT/AD[0:0] + FMT/AD[0:1] >= 10 & FMT/AD[0:1] >=3)'
echo -e "Filtering Expression: $filter_expression"

for vcf in ../../data/vcf_filtermutectcalls_obf/*/*.vcf; do
    
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