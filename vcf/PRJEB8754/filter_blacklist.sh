#!/bin/bash

set -euo pipefail

root_outdir="vcf_filtered_pass-orient-pos-sb-ad-blacklist"
mkdir -p $root_outdir

blacklist_path="../../data/blacklists/master_blacklist.bed.gz"
echo -e "Filtering using blacklist: $blacklist_path"

for vcf in vcf_filtered_pass-orient-pos-sb-ad_dup-unmarked/*/*.vcf; do
    
    filename=$(basename $vcf)
    sample_name=${filename%%.*}

    echo $sample_name
    
    outdir="${root_outdir}/${sample_name}"
    mkdir -p $outdir

    echo -e "\nFiltering $filename"
    bcftools view -T ^"$blacklist_path" $vcf -o "${outdir}/${sample_name}.vcf"
    echo "Filtered VCF saved to ${outdir}/${sample_name}.vcf"

    echo "Indexing ${outdir}/${sample_name}.vcf"
    gatk IndexFeatureFile -I "${outdir}/${sample_name}.vcf"

done

echo -e "\nFinished."