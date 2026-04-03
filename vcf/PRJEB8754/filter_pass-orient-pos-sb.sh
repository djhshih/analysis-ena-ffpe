#!/bin/bash

set -euo pipefail

root_indir="dup-unmarked_raw_filtermutectcalls_obmm_unfiltered"
root_outdir="dup-unmarked_filtered_pass-orient-pos-sb"
mkdir -p $root_outdir

filter_expression='(FILTER="PASS" | FILTER="strand_bias" | FILTER="orientation" | FILTER="position" | FILTER="strand_bias;orientation" | FILTER="strand_bias;position" | FILTER="orientation;position" | FILTER="strand_bias;orientation;position")'
echo -e "Filtering Expression: $filter_expression"

for vcf in $root_indir/*/*.vcf; do
    
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
