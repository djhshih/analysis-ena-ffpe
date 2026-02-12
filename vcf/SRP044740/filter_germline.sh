#!/usr/bin/env bash

set -euo pipefail

target_indir="../../germline-filter/SRP044740/filtered_pass-orientation-dp20-blacklist-macni"
outdir_root=filtered_pass-orientation-dp20-blacklist-macni

for vcf in filtered_pass-orientation-dp20-blacklist/*/*.vcf; do

    echo Filtering: $vcf

    sample_name=$(basename $(dirname $vcf))
    target_list=$target_indir/$sample_name/$sample_name.tsv

    if [[ -f $target_list ]]; then

        outdir=$outdir_root/$sample_name
        mkdir -p $outdir

        bcftools view -T $target_list $vcf -o $outdir/$sample_name.vcf

    else
        echo $target_list does not exist. 
        echo Likely reason is that all variants were filtered out

    fi

done

echo Outputs saved to: $outdir_root