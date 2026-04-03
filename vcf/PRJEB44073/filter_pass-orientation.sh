#!/bin/bash

set -euo pipefail

indir_root="raw_filtermutectcalls-obmm"
outdir_root="filtered_pass-orientation"
mkdir -p $outdir_root

filter_expression='(FILTER="PASS" | FILTER="orientation")'
echo -e "Filtering Expression: $filter_expression"

i=1
for vcf in $indir_root/*/*.vcf; do
    
    filename=$(basename $vcf)
    sample_name=${filename%%.*}

    echo $i. Filtering Sample: $sample_name
    
    outdir="${outdir_root}/${sample_name}"
    mkdir -p $outdir
    outpath="${outdir}/${sample_name}.vcf.gz"

    bcftools view -i "$filter_expression" $vcf -o "$outpath"
    echo -e "\t- Filtered VCF saved to $outpath"

    bcftools index -t "$outpath"
    echo -e "\t- Indexed $outpath"

    ((i++))
done

echo -e "\nFinished."
