#!/usr/bin/env bash

set -euo pipefail

outdir_root=dup-unmarked_filtered_pass-orient-pos-sb-vaf-dp-macni
target_indir="../../germline-filter/PRJEB8754/dup-unmarked_filtered_pass-orient-pos-sb-vaf-dp-macni"


for vcf in dup-unmarked_filtered_pass-orient-pos-sb-ad-blacklist/*/*.vcf.gz; do

    echo $vcf

    sample_name=$(basename $(dirname $vcf))
    target_list=$target_indir/$sample_name/$sample_name.tsv

    if [[ -f $target_list ]]; then

        outdir=$outdir_root/$sample_name
        mkdir -p $outdir

        bcftools view -T $target_list $vcf -Oz -o $outdir/$sample_name.vcf.gz
        bcftools index -t $outdir/$sample_name.vcf.gz

    else
        echo $target_list does not exist. 
        echo Find out why

    fi

done